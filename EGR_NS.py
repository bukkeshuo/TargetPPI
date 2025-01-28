import torch
import os
import time
import sys
import numpy as np

import torch as t
from torch import nn
from torch.autograd import Variable
import torch.nn.functional as F

sys.path.append("../")
from config import DefaultConfig,parameter

configs = DefaultConfig()
window_size=parameter['win_size']
# print("window_size",window_size)

def init_weights(m):#更适用于使用ReLU及其变种作为激活函数的网络结构。
    # 因为ReLU在正区间内梯度恒为1，He初始化有助于保持这一特性，避免在训练过程中梯度消失或爆炸。
    if isinstance(m, nn.Linear):
        torch.nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
        if m.bias is not None:
            m.bias.data.fill_(0.01)
    elif isinstance(m, nn.Sequential):
        for layer in m:
            init_weights(layer)  # 递归地对每层应用初始化

def initialize_weight(x):#适用于使用sigmoid或tanh等S型激活函数的网络结构。这些激活函数在输入为0时导数最大，
    # 因此使用Xavier初始化可以确保网络中各层的激活输入在训练初期保持在有效的导数范围内。
    if isinstance(x, nn.Linear):
        nn.init.xavier_uniform_(x.weight)
        if x.bias is not None:
            nn.init.constant_(x.bias, 0)
    elif isinstance(x, nn.Sequential):
        for layer in x:
            initialize_weight(layer)  # 递归地对每层应用初始化


class egret_ppi(nn.Module):
    def __init__(self, ratio=None):
        super(egret_ppi, self).__init__()
        global configs


        ####包含两个层的序列：一个线性层和一个 Sigmoid 激活函数。这个序列用于映射模型的输出到二进制分类的概率。
        self.outLayer = nn.Sequential(
            nn.Linear(32+64+32,1),
            nn.Sigmoid())

        initialize_weight(self.outLayer)

        self.t5con1d = nn.Sequential(
            nn.Conv1d(in_channels=1024,
                      out_channels=64,
                      kernel_size=21, stride=1,
                      padding=21 // 2, dilation=1, groups=1,
                      bias=True, padding_mode='zeros'),
            nn.LeakyReLU(negative_slope=.01),
            nn.BatchNorm1d(num_features=64,
                           eps=1e-05, momentum=0.1, affine=True, track_running_stats=True),
            nn.Dropout(.5)  # it was .2
        )

        self.layer_norm = nn.LayerNorm(normalized_shape=33,eps=1e-6)
        self.gat_norm = nn.LayerNorm(normalized_shape=36,eps=1e-6)
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        ###卷积神经网络（CNN）编码器。这个编码器包括了一个卷积层、LeakyReLU 激活函数、批归一化和 dropout。
        # 卷积层的输入通道数为 1024，输出通道数为 32，卷积核大小为 7，使用 LeakyReLU 作为激活函数，并在之后应用批归一化和 dropout。
        self.relu = nn.ReLU()

        self.bidirectional_lstm_feature = nn.LSTM(input_size=33, hidden_size=16, num_layers=1, batch_first=True,
                                             bidirectional=True)

        self.dense_feature = nn.Linear((window_size*2+1) * 32,32)  # 64*2 from bidirectional LSTM

        import EAGA_NS as egret
        #model/multiegrt.py
        config_dict = egret.config_dict

        self.gat_layer = egret.GAT(
            # in_dim=36+256+256,
            in_dim=32 + 64,
            hidden_dim=64,
            out_dim=32,
            edge_dim=6,
            num_heads=1,
            use_bias=False,
            merge1='cat',  # 第一层的 merge 参数
            merge2='mean',  # 第二层的 merge 参数
            config_dict=config_dict)



    def forward(self, feature,t5residue,graph_batch):#protbert_feature 是输入的蛋白质特征。graph_batch 是输入的图数据。

        t5residue = t5residue.permute(0, 2, 1)
        t5residue = self.t5con1d(t5residue)
        t5residue = t5residue.permute(0, 2, 1)

        batch_size,seq_lenth,local_dim,feature_dim = feature.shape
        y=z=self.layer_norm(feature)
        y=y.view(batch_size*seq_lenth,local_dim,feature_dim)
        y, _ = self.bidirectional_lstm_feature(y)
        y = y.reshape(batch_size, seq_lenth, -1)
        y = F.relu(self.dense_feature(y))

        z = torch.cat((y, t5residue), dim=2)

        features2, head_attn_scores = self.gat_layer(graph_batch, z.view([batch_size*seq_lenth, 32 + 64]))
        features2 = features2.view([batch_size, seq_lenth, 32])
        features = t.cat((features2, z), 2)
        features = self.outLayer(features)

        return features, head_attn_scores


