#coding=UTF-8

from config import parameter
import torch
from EGR_NS import egret_ppi

from DataProcess import dataSet,graph_collate
from torch.utils.data import DataLoader
from sklearn.metrics import roc_curve, precision_recall_curve,auc

import os

os.environ["CUDA_VISIBLE_DEVICES"] = "0"
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device)

import numpy as np
import math
def CalculateEvaluationMetrics(y_true, y_pred):
    TP = float(0)
    FP = float(0)
    TN = float(0)
    FN = float(0)
    for i, j in zip(y_true, y_pred):
        if (i == 1 and j == 1):
            TP += 1
        elif (i == 0 and j == 1):
            FP += 1
        elif (i == 0 and j == 0):
            TN += 1
        elif (i == 1 and j == 0):
            FN += 1
    print("TP: ", TP)
    print("FP: ", FP)
    print("TN: ", TN)
    print("FN: ", FN)
    sensitivity = TP / (TP + FN)
    specificity = TN / (TN + FP)
    recall = TP / (TP + FN)
    precision = TP / (TP + FP)
    MCC = (TP * TN - FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    F1 = 2 * (precision * recall) / (precision + recall)
    accuracy = (TP+TN)/(TP+TN+FP+FN)
    print("accuracy:",accuracy)
    print("sensitivity: ", sensitivity)
    print("Specificity: ", specificity)
    # same as sensitivity
    print("Recall: ", recall)
    print("Precision: ", precision)
    print("MCC: ", MCC)
    print("F1: ", F1)
    return TP, FP, TN, FN, sensitivity, specificity, recall, precision, MCC, F1, accuracy

def Compute_data(models):

    result_out = r'存放结果的txt'
    result_lab = r'存放lable的txt'

    file_out = open(result_out, 'w')
    file_lab = open(result_lab, 'w')

    for batch in validation_loader:
        feature, localt5residue, labels, pro_lengths, graph_batch,pro_name = batch

        feature = torch.tensor([item.cpu().detach().numpy() for item in feature]).cuda()
        labels = torch.tensor([item.cpu().detach().numpy() for item in labels]).cuda()
        localt5residue = torch.tensor([item.cpu().detach().numpy() for item in localt5residue]).cuda()

        pro_name_str = pro_name[0]

        graph_batch = graph_batch.to(device)
        out_num=0
        graph_batch.edata['ex'] = graph_batch.edata['ex'].to(device).float()
        for model_index in range(len(models)):
            model = models[model_index]
            output, _ = model(feature, localt5residue, graph_batch)
            out_num+=output
        out=out_num/len(models)

        out = out.view(-1)
        labels=labels.view(-1)

        for out_i, label_value in zip(out, labels):
            file_out.write(f"{out_i.item()}\n")
            file_lab.write(f"{label_value.item()}\n")
    file_out.close()
    file_lab.close()



    # 读取预测结果文件和真实标签文件
    with open(result_out, 'r') as file_out, open(result_lab, 'r') as file_lab:
        content_result = file_out.readlines()
        content_label = file_lab.readlines()
    # 将内容转换为浮点数列表
    pred = [float(item.strip()) for item in content_result]
    # print(pred)
    truth = [float(item.strip()) for item in content_label]

    # 计算阈值
    sorted_pred = np.sort(pred)
    sorted_pred_descending = np.flip(sorted_pred)
    num_of_1 = np.count_nonzero(truth)
    print(num_of_1)
    threshold = sorted_pred_descending.item(num_of_1 - 1)
    print(threshold)
    # pred_binary_sum = sum(list(np.where(pred > threshold, 1, 0)))
    pred_binary_sum = sum([1 if x > threshold else 0 for x in pred])
    print(pred_binary_sum)
    pred_binary = []
    flag = 0
    for item in pred:
        if item == threshold:
            if flag < num_of_1-pred_binary_sum:
                pred_binary.append(1)
                flag+=1
            else:
                pred_binary.append(0)
        elif item > threshold:
            pred_binary.append(1)
        else:
            pred_binary.append(0)

    TP, FP, TN, FN, sensitivity, specificity, recall, precision, MCC, F1_score, accuracy = CalculateEvaluationMetrics(truth, pred_binary)
    # PrintToCSV(csvPre + args_prefix, au_roc, aupr, TP, FP, TN, FN, sensitivity, specificity, recall, precision, MCC,
    #                threshold, F1_score, accuracy)


    from sklearn.metrics import roc_auc_score, average_precision_score
    auc_score = roc_auc_score(truth, pred)

    # 计算AUPR
    aupr_score = average_precision_score(truth, pred)
    print("AUC:", auc_score)
    print("AUPR:", aupr_score)



num_list = ['164','186','72','448','355','70','60']

'+num+'
for num in num_list:
    print(num)
    feature_path="you path"

    protein_path = feature_path+'/dset_'+num
    pssm_path = feature_path+'/pssm_'+num
    hydropathy_path = feature_path+'/hydropathy_'+num
    Physical_properties_path = feature_path+'/Physical_properties_'+num
    Physiochemical_characteristics_path = feature_path+'/Physiochemical_characteristics_'+num
    pKx_path = feature_path+'/pKx_'+num
    position_path = feature_path+'/position_'+num

    t5_path=feature_path+num+'/t5'#t5路径

    protein_list_file = feature_path+num+'.txt'#测试数据集 蛋白名

    dist_matrix_path=feature_path+'/dist_mat'
    angle_matrix_path=feature_path+'/angle_mat'

    Get_Data_Validation=dataSet(pssm_path, hydropathy_path, Physical_properties_path,
                 Physiochemical_characteristics_path, pKx_path, position_path,
                 protein_list_file, protein_path, dist_matrix_path, angle_matrix_path,t5_path)

    validation_loader = DataLoader(dataset=Get_Data_Validation, batch_size=1, shuffle=False,collate_fn=graph_collate)

    models = []
    model_path = r'model_path'
    for i in range(9):
        model = egret_ppi(ratio=None).to(device)
        model_file = f'{model_path}{i}model1'
        try:
            model.load_state_dict(torch.load(model_file, map_location='cpu'), strict=True)  # 加载模型状态
        except FileNotFoundError:
            print(f"Model file {model_file} not found.")
            continue
        model.eval()
        models.append(model)  #
    result = Compute_data(models)
