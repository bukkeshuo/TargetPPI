import os
import time
import pickle
import torch as t
import numpy as np
from torch.utils.data import Dataset, DataLoader
import gzip
from time import time
from config import DefaultConfig, parameter
import torch
import dgl

LOCAL_NEGHBOR_SIZE = parameter['win_size']
class ProteinProcessor:
    def __init__(self, features):
        if isinstance(features, torch.Tensor):
            self.features = features.numpy()
        else:
            self.features = features
        self.proteinlen = self.features.shape[0]

    def read_features_2d(self, neighbor_list):
        selected_feature = []
        for neighbor in neighbor_list:
            if neighbor != 'Zpad':
                try:
                    selected_feature.append(self.features[neighbor])
                except IndexError:
                    exit(1)
            else:
                selected_feature.append(np.zeros(self.features.shape[1]))
        return np.array(selected_feature)

    def build_2d_windows(self, window_size):
        LOCAL_NEIGHBOR_SIZE = window_size
        all_features_2d = []
        for aa_index in range(self.proteinlen):
            add_to_left = max(0, LOCAL_NEIGHBOR_SIZE - aa_index)
            add_to_right = max(0, aa_index + LOCAL_NEIGHBOR_SIZE + 1 - self.proteinlen)

            neighbor_list = [i for i in range(aa_index - LOCAL_NEIGHBOR_SIZE + add_to_left,
                                              aa_index + LOCAL_NEIGHBOR_SIZE + 1 - add_to_right)]

            cnt = 0
            lr_flag = add_to_right - add_to_left
            while cnt < add_to_right + add_to_left:
                if lr_flag <= 0:
                    neighbor_list.insert(0, 'Zpad')
                    lr_flag += 1
                else:
                    neighbor_list.append('Zpad')
                    lr_flag -= 1
                cnt += 1

            features_2d = self.read_features_2d(neighbor_list)
            all_features_2d.append(features_2d)

        return np.array(all_features_2d)

class dataSet(Dataset):
    def __init__(self, pssm_path, hydropathy_path, Physical_properties_path,
                 Physiochemical_characteristics_path, pKx_path, position_path,
                 protein_list_file, protein_path, dist_matrix_path, angle_matrix_path,t5_path):
        super(dataSet, self).__init__()


        self.edge_feat_mean = [40.02347537780368, 1.640517914996531, 39.73830284206693, 1.5588000279344325, 1.5438419537544195, 1.5698044805435276]
        self.edge_feat_std = [26.48078921810372, 0.48689876318729164, 26.496414431035085, 0.6903461210048146, 0.9065970028521069, 0.6831223214174721]


        self.pssm_path = pssm_path
        self.protein_path = protein_path
        self.hydropathy_path = hydropathy_path
        self.Physical_properties_path = Physical_properties_path
        self.Physiochemical_characteristics_path = Physiochemical_characteristics_path
        self.pKx_path = pKx_path
        self.position_path = position_path
        self.dist_matrix_path = dist_matrix_path
        self.angle_matrix_path = angle_matrix_path
        self.t5_path = t5_path



        with open(protein_list_file, "r") as f:
            protein_list = f.readlines()
            self.protein_list = [x.strip() for x in protein_list]


        self.neighbourhood_size = 21

        self.protein_list_len = len(self.protein_list)

    def __getitem__(self, index):
        pro_name = self.protein_list[index]
        id_idx = index

        feature,t5residue,label, pro_length = self.GetResidues(pro_name)

        G = self.generate_graph(pro_length,pro_name)

        return feature,t5residue,label, pro_length, G,pro_name

    def __len__(self):
        return self.protein_list_len

    def generate_graph(self,pro_length,protein_name):
        G = dgl.DGLGraph()
        G.add_nodes(pro_length)

        dist_mat_path = os.path.join(self.dist_matrix_path, f"{protein_name}_dist.npy")
        angle_mat_path = os.path.join(self.angle_matrix_path, f"{protein_name}_angle.npy")

        dist_mat = np.load(dist_mat_path)
        angle_mat = np.load(angle_mat_path)

        neighborhood_indices = dist_mat[:, :, 0] .argsort()[:, 1:self.neighbourhood_size]

        edge_feat = np.array([
            dist_mat[:, :, 0],
            dist_mat[:, :, 2],
            dist_mat[:, :, 3],
            angle_mat[:, :, 0],
            angle_mat[:, :, 1],
            angle_mat[:, :, 2]
        ])


        edge_feat = np.transpose(edge_feat, (1, 2, 0))
        edge_feat = (edge_feat - self.edge_feat_mean) / self.edge_feat_std

        self.add_edges_custom(G, neighborhood_indices, edge_feat)

        return G

    def add_edges_custom(self, G, neighborhood_indices, edge_features):
        size = neighborhood_indices.shape[0]
        neighborhood_indices = neighborhood_indices.tolist()
        src = []
        dst = []
        temp_edge_features = []

        for center in range(size):
            src += neighborhood_indices[center]
            dst += [center] * (self.neighbourhood_size - 1)
            for nbr in neighborhood_indices[center]:
                temp_edge_features += [np.abs(edge_features[center, nbr])]

        if len(src) != len(dst):
            print('source and destination array should have been of the same length: src and dst:', len(src), len(dst))
            raise Exception
        G.add_edges(src, dst)
        G.edata['ex'] = torch.tensor(temp_edge_features)

    def GetResidues(self, pro_name):


        pssm_name_path = os.path.join(self.pssm_path, pro_name + '.opssm')
        pssm = np.genfromtxt(pssm_name_path, skip_footer=0, skip_header=0)[:, :20]
        pro_length, psa_dim = pssm.shape

        hydropathy_name_path = os.path.join(self.hydropathy_path, pro_name + '.txt')
        hydropathy = np.genfromtxt(hydropathy_name_path, skip_footer=0, skip_header=0).reshape(pro_length, 1)

        Physical_properties_name_path = os.path.join(self.Physical_properties_path, pro_name + '.txt')
        Physical_properties = np.genfromtxt(Physical_properties_name_path, skip_footer=0, skip_header=0)

        Physiochemical_characteristics_name_path = os.path.join(self.Physiochemical_characteristics_path,
                                                                pro_name + '.txt')
        Physiochemical_characteristics = np.genfromtxt(Physiochemical_characteristics_name_path, skip_footer=0,
                                                       skip_header=0)
        pro_length, Physiochemical_characteristics_dim = Physiochemical_characteristics.shape

        pKx_name_path = os.path.join(self.pKx_path, pro_name + '.txt')
        pKx = np.genfromtxt(pKx_name_path, skip_footer=0, skip_header=0).reshape(pro_length, 1)

        position_name_path = os.path.join(self.position_path, pro_name + '.txt')
        position = np.genfromtxt(position_name_path, skip_footer=0, skip_header=0).reshape(pro_length, 1)
        # print("protenlenth",pro_length)

        fea_pssm = torch.FloatTensor(pssm)
        fea_hydropathy = torch.FloatTensor(hydropathy)
        fea_Physical_properties = torch.FloatTensor(Physical_properties)
        fea_Physiochemical_characteristics = torch.FloatTensor(Physiochemical_characteristics)
        fea_pKx = torch.FloatTensor(pKx)
        fea_position = torch.FloatTensor(position)
        new_feature_nopsa = torch.cat((fea_pssm, fea_hydropathy, fea_Physical_properties,
                                 fea_Physiochemical_characteristics, fea_pKx, fea_position), 1)#33


        new_feature_nopsa= new_feature_nopsa

        processor = ProteinProcessor(new_feature_nopsa)
        features_3d = processor.build_2d_windows(LOCAL_NEGHBOR_SIZE)
        features_3d = torch.FloatTensor(features_3d)

        t5_path = os.path.join(self.t5_path, f"{pro_name}.embd")
        t5residue = np.genfromtxt(t5_path, skip_footer=0, skip_header=0)
        t5residue = torch.FloatTensor(t5residue)

        file = open(os.path.join(self.protein_path, pro_name))
        label = file.readline()
        label = file.readline()
        label = file.readline()
        file.close()
        label = list(label.strip())
        label = list(map(int, label))
        label = np.array(label)
        label_length = len(label)
        label = label.reshape(len(label), 1)

        return features_3d,t5residue,label, label_length


def graph_collate(data):
    features,t5residue,labels, pro_lengths, graphs,pro_name = zip(*data)

    features = [torch.tensor(feature, dtype=torch.float32) if isinstance(feature, np.ndarray) else feature for feature
                in features]
    t5residues = [torch.tensor(residue, dtype=torch.float32) if isinstance(residue, np.ndarray) else residue for
                       residue in t5residue]

    labels = [torch.tensor(label, dtype=torch.float32) if isinstance(label, np.ndarray) else label for label in labels]
    pro_lengths = torch.tensor(pro_lengths, dtype=torch.int64)  # 根据实际情况选择适当的数据类型
    t5residues = torch.stack(t5residues)
    graphs = dgl.batch(graphs)
    pro_name = list(pro_name)

    return features,t5residues, labels, pro_lengths, graphs,pro_name
