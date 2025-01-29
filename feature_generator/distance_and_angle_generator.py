# coding: utf-8

import numpy as np
from time import time
import gzip
import warnings
import pickle
from Bio.PDB import *
import os


warnings.filterwarnings("ignore")


parser = PDBParser()
THIRD_ATOM = 'N'  # 'O'


def residue_distance(res1, res2):
    distance = []
    cnt = 0
    for atom1 in res1:
        for atom2 in res2:
            distance += [abs(atom1 - atom2)]
            cnt += 1
    distance = np.array(distance)
    dist_mean = distance.mean()
    dist_std = distance.std()
    if 'CB' in res1 and 'CB' in res2:
        dist_cb = abs(res1['CB'] - res2['CB'])
    else:
        dist_cb = dist_mean
    if 'CA' in res1 and 'CA' in res2:
        dist_ca = abs(res1['CA'] - res2['CA'])
    else:
        dist_ca = dist_mean
    return dist_mean, dist_std, dist_cb, dist_ca


def residue_relative_angle(res1, res2):
    if 'CA' in res1 and THIRD_ATOM in res1 and 'C' in res1:
        v1 = res1['CA'].get_vector().get_array()
        v2 = res1[THIRD_ATOM].get_vector().get_array()
        v3 = res1['C'].get_vector().get_array()
        a1 = np.cross((v2 - v1), (v3 - v2))
    else:
        k = list(res1)
        if len(k) == 1:
            a1 = k[0].get_vector().get_array()
        else:
            raise
    a1 = a1 / np.linalg.norm(a1)

    if 'CA' in res2 and THIRD_ATOM in res2 and 'C' in res2:
        v1 = res2['CA'].get_vector().get_array()
        v2 = res2[THIRD_ATOM].get_vector().get_array()
        v3 = res2['C'].get_vector().get_array()
        a2 = np.cross((v2 - v1), (v3 - v2))
    else:
        k = list(res2)
        if len(k) == 1:
            a2 = k[0].get_vector().get_array()
        else:
            raise
    a2 = a2 / np.linalg.norm(a2)
    a = np.arccos(np.clip(np.dot(a1, a2), -1.0, 1.0))

    if 'CA' in res1 and 'CB' in res1 and 'CB' in res2:
        v1 = res1['CA'].get_vector().get_array()
        v2 = res1['CB'].get_vector().get_array()
        v3 = res2['CB'].get_vector().get_array()
        w1 = np.cross((v2 - v1), (v3 - v2))
    elif 'CB' not in res1 or 'CB' not in res2:
        k = list(res1)
        w1 = k[0].get_vector().get_array()
    else:
        k1 = list(res1)
        k2 = list(res2)
        if len(k1) == 1:
            w1 = k1[0].get_vector().get_array()
        elif len(k2) == 1 and len(k1) != 1:
            v4 = k2[0].get_vector().get_array()
            w1 = np.cross((v2 - v1), (v4 - v2))
        else:
            raise
    w1 = w1 / np.linalg.norm(w1)

    if 'CA' in res2 and 'CB' in res2 and 'CB' in res1:
        v1 = res2['CA'].get_vector().get_array()
        v2 = res2['CB'].get_vector().get_array()
        v3 = res1['CB'].get_vector().get_array()
        w2 = np.cross((v2 - v1), (v3 - v2))
    elif 'CB' not in res1 or 'CB' not in res2:
        k = list(res2)
        w2 = k[0].get_vector().get_array()
    else:
        k1 = list(res1)
        k2 = list(res2)
        if len(k1) == 1 and len(k2) != 1:
            v4 = k1[0].get_vector().get_array()
            w2 = np.cross((v2 - v1), (v4 - v2))
        elif len(k2) == 1:
            w2 = k2[0].get_vector().get_array()
    w2 = w2 / np.linalg.norm(w2)
    w = np.arccos(np.clip(np.dot(w1, w2), -1.0, 1.0))

    if 'CA' in res1 and THIRD_ATOM in res1 and 'CB' in res1:
        v1 = res1['CA'].get_vector().get_array()
        v2 = res1[THIRD_ATOM].get_vector().get_array()
        v3 = res1['CB'].get_vector().get_array()
        theta1 = np.cross((v2 - v1), (v3 - v2))
    elif 'CB' not in res1:
        k = list(res1)
        theta1 = k[0].get_vector().get_array()
    else:
        k1 = list(res1)
        if len(k1) == 1:
            theta1 = k1[0].get_vector().get_array()
        else:
            raise
    theta1 = theta1 / np.linalg.norm(theta1)

    if 'CA' in res2 and 'CB' in res2 and 'CB' in res1:
        v1 = res2['CA'].get_vector().get_array()
        v2 = res2['CB'].get_vector().get_array()
        v3 = res1['CB'].get_vector().get_array()
        theta2 = np.cross((v2 - v1), (v3 - v2))
    elif 'CB' not in res1 or 'CB' not in res2:
        k = list(res2)
        theta2 = k[0].get_vector().get_array()
    else:
        k1 = list(res1)
        k2 = list(res2)
        if len(k1) == 1 and len(k2) != 1:
            v4 = k1[0].get_vector().get_array()
            theta2 = np.cross((v2 - v1), (v4 - v2))
        elif len(k2) == 1:
            theta2 = k2[0].get_vector().get_array()
    theta2 = theta2 / np.linalg.norm(theta2)
    theta = np.arccos(np.clip(np.dot(theta1, theta2), -1.0, 1.0))
    return a, w, theta


def get_dist_and_angle_matrix(residues):
    size = len(residues)
    dist_mat = np.zeros([size, size, 4])
    angle_mat = np.zeros([size, size, 3])
    for i in range(size):
        for j in range(i + 1, size):
            dist_mean, dist_std, dist_ca, dist_c = residue_distance(residues[i], residues[j])
            a, w, theta = residue_relative_angle(residues[i], residues[j])

            dist_mat[i, j, 0] = dist_mean
            dist_mat[i, j, 1] = dist_std
            dist_mat[i, j, 2] = dist_ca
            dist_mat[i, j, 3] = dist_c

            dist_mat[j, i, 0] = dist_mean
            dist_mat[j, i, 1] = dist_std
            dist_mat[j, i, 2] = dist_ca
            dist_mat[j, i, 3] = dist_c

            angle_mat[i, j, 0] = a
            angle_mat[i, j, 1] = w
            angle_mat[i, j, 2] = theta

            angle_mat[j, i, 0] = a
            angle_mat[j, i, 1] = w
            angle_mat[j, i, 2] = theta
    return dist_mat, angle_mat


from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqUtils

protein_letters_3to1 = SeqUtils.IUPACData.protein_letters_3to1
ppb = PPBuilder()


def generate_distance_and_angle_matrix(root_dir):
    t0 = time()
    with open(root_dir, 'r') as f:
        protein_list = f.readlines()
        protein_list = [x.strip() for x in protein_list]



    save_path1 = r"your path/dist_mat"
    if not os.path.exists(save_path1):
        os.makedirs(save_path1)
    save_path2 = r"your pdb path/angle_mat"
    if not os.path.exists(save_path2):
        os.makedirs(save_path2)

    # for idx, protein_name in enumerate(protein_list):
    for idx, protein_name in enumerate(protein_list):
        dist_mat_path = os.path.join(save_path1, f"{protein_name}_dist.npy")
        angle_mat_path = os.path.join(save_path2, f"{protein_name}_angle.npy")

        try:
            path = r'your pdb path/' + protein_name + '.pdb'

            # structure=parser.get_structure(path)
            structure = parser.get_structure(protein_name, path)
            print(f"Structure for {protein_name} parsed successfully.")
            model = structure[0]  # every structure object has only one model object in DeepPPISP's dataset
            pep_pdb = ''
            residues = []

            for residue in model.get_residues():
                print(residue)
                residues += [residue]
            print(f"Total number of residues: {len(residues)}")
            peptides = ppb.build_peptides(model)
            pep_pdb += ''.join([str(pep.get_sequence()) for pep in peptides])

            dist_mat, angle_mat = get_dist_and_angle_matrix(residues)
            # 保存每个蛋白质的距离矩阵和角度矩阵文件
            np.save(dist_mat_path, dist_mat)
            np.save(angle_mat_path, angle_mat)

            # print(f"{idx  + 1}/{len(protein_list)}: {protein_name} processed.")
        except Exception as e:
            print(f"Error processing {protein_name}: {e}")  # 打印详细错误信息
            print("error")# print(f"{idx + 1}/{len(protein_list)}: {protein_name} failed with error: {e}")

    print(f"Total time: {time() - t0} seconds.")


generate_distance_and_angle_matrix("/home/dm/EGRET-Transform/inputs/dset70.txt")#the file include protein name which you want generate

