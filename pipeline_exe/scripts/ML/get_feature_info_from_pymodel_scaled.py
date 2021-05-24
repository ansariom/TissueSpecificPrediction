'''
Created on Oct 12, 2017

@author: mitra
'''

import sys, getopt
import pickle
import numpy as np
import pandas as pd
import csv
from sklearn import preprocessing
from sklearn import metrics
from sklearn.metrics import average_precision_score
from numpy import int32



'================================================'
'   run saved model on given set of features  '
'================================================'
def run_model(x_test, x_test_tss_id, outfile, model_file_name, class_labels = None):
    trained_model = pickle.load(open(model_file_name, 'rb'))
    result = trained_model.predict_proba(x_test)
    res_records = np.rec.fromarrays((x_test_tss_id, result[:,0].ravel(), result[:,1].ravel()), names= ('tss_name', 'prob0', 'prob1'))
    df = pd.DataFrame.from_records(res_records)
    df.to_csv(outfile, sep = "\t", columns = ('tss_name', 'prob0', 'prob1'), index = False, header = True)
    
    if class_labels != None:
    	fpr, tpr, thresholds = metrics.roc_curve(class_labels, result[:,1], pos_label=1)
    	auroc = metrics.auc(fpr, tpr)
    	average_precision = average_precision_score(class_labels, result[:,1])
    	print(auroc, average_precision)

if __name__ == '__main__':
    args = sys.argv[1:]
    model_file = args[0] # saved model
    train_set = args[1]  # unscaled train array  > traind.npy
    all_features_csv = args[2]    # features.csv without class labels
    outfile = args[3]    # outfile for probability of tissue-spec for given list of features
    scaled_outfile = args[4]  # output the scaled array for given all_features_csv
    scaling = args[5]    #1=no_scale, 2=mean and sd 3=0-1
   
    #print all_features_csv 
    # read input features, 1st column is tss_id and last column is class label
    x_test = []
    first_line = True
    with open(all_features_csv, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            if (first_line):
                features_header = np.asarray(row)
                first_line = False
                continue    
            x_test.append(np.asarray(row))
    x_test = np.asanyarray(x_test)
    #x_test = pd.read_csv(all_features_csv)

    # 1- compute all_features_scaled using saved train_set 
    print(x_test.shape)
    x_test_features = np.asanyarray(x_test[:,1:], np.float32)
    print(x_test_features.shape)
    x_test_tss_ids = np.asarray(x_test[:,0])
    print(x_test_tss_ids.shape)
    x_train = np.load(train_set)
    scaler = preprocessing.StandardScaler(with_mean=True, with_std=True).fit(x_train)
    #scaler = preprocessing.MinMaxScaler().fit(x_train) 
    x_test_scaled = scaler.transform(x_test_features)
    
    run_model(x_test_scaled, x_test_tss_ids, outfile, model_file)
    #run_model(x_test_scaled, x_test_tss_ids, outfile, model_file, np.asanyarray(x_test[:,-1], int32))
    # Save scaled dataset for all TSS
    features_header = np.concatenate((["tss_name"], features_header), axis = 0)
    x_test_tss_ids = np.asarray(x_test[:,0]).reshape(1,-1)
    scaled_data = np.concatenate((x_test_tss_ids.T, x_test_scaled), axis = 1)
    df = pd.DataFrame(scaled_data, columns = np.asarray(features_header))
    df.to_csv(scaled_outfile, sep = "\t", index = False, header = True)
    
    
    
    
