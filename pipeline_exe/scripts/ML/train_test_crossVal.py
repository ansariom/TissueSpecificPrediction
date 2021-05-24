import sys, getopt
import numpy as np
import csv
from sklearn import preprocessing
from numpy import float32, int32
from sklearn.preprocessing.data import StandardScaler
from sklearn import metrics
from sklearn import linear_model
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import average_precision_score
from os.path import os
import pickle
from __builtin__ import str

max_param = 0.01
all_features = []
test_all_features = []
train_all_labels = []
test_all_labels = []
train_all_features = []
features_header = []
batch_index = 0
learning_rate = 0.01
training_epochs = 200
batch_size = 100
display_step = 50
weight_decays = np.arange(0.0001, max_param ,0.002) 
fold = 0
fold_results = []
outdir = "-"

'Split data into train and helout test'
def split_data_train_test(seed):
    global train_all_features, test_all_features, train_all_labels, test_all_labels
    print(seed)
    nfeatures = all_features.shape[0]
    ntest = int(0.2 * nfeatures)
    np.random.seed(seed)
    
    test_feature_ids = np.random.choice(nfeatures, ntest, replace = False)
    test_all = all_features[test_feature_ids]
    
    all_idx = np.arange(nfeatures)
    train_ids = np.setdiff1d(all_idx, test_feature_ids)
    train_all = all_features[train_ids]
    
    train_all_labels = np.asanyarray(train_all[:,-1], int32).reshape(len(train_all), 1)
    train_all_features = np.asanyarray(train_all[:,1:-1], float32)
    
    test_all_labels = np.asanyarray(test_all[:,-1], int32).reshape(len(test_all), 1)
    test_all_features = np.asanyarray(test_all[:,1:-1], float32)

    train_tss_list = outdir + "/" + argsDict["model_type"] + "_train_tss_ids.txt"
    test_tss_list = outdir + "/" + argsDict["model_type"] + "_test_tss_ids.txt"
    np.savetxt(train_tss_list, train_all[:,0], delimiter = "\t", fmt="%s")
    np.savetxt(test_tss_list, test_all[:,0], delimiter = "\t" , fmt = "%s")
    
def split_train_cross_val(nfold):
    'split train into a pie with 5 section, 4 train 1 validation'
    param_file = outdir + "/" + argsDict["model_type"] + "_crossval_params_performace.txt"
    f = open(param_file, "wt")
    
    fold_results = []
    ntrain = train_all_features.shape[0]
    pie_size = ntrain/nfold
    all_idx = np.arange(ntrain)
    
    for i in xrange(nfold):
        fold = i
        start_idx = i * pie_size
        end_idx = start_idx + pie_size
        if (i == nfold - 1):
            end_idx = ntrain
        test_fold_features = train_all_features[start_idx:end_idx]
        test_fold_labels = train_all_labels[start_idx:end_idx]
        test_fold_labels = test_fold_labels.reshape(len(test_fold_labels), 1)
        
        train_idxs = np.setdiff1d(all_idx, np.arange(start_idx, end_idx))
        train_fold_features = train_all_features[train_idxs]
        train_fold_labels = train_all_labels[train_idxs]
        train_fold_labels = train_fold_labels.reshape(len(train_fold_labels), 1)
        
        # normalize the data
#        scaler = preprocessing.MinMaxScaler().fit(train_fold_features)
        scaler = preprocessing.StandardScaler(with_mean=True, with_std=True).fit(train_fold_features)
        #scaler = preprocessing.RobustScaler().fit(train_fold_features)
        scaled_train_features = scaler.transform(train_fold_features)
        scaled_test_features = scaler.transform(test_fold_features)
        for param in iter(weight_decays):
            auroc, auprc = linear_model_simple(scaled_train_features, train_fold_labels, scaled_test_features, test_fold_labels, param)
            fold_results.append(np.array([fold, param, auroc, auprc]))
            f.write(str(fold) + "\t" + str(param) + "\t" + str(auroc) + "\t" + str(auprc) + "\n")
#             print(fold, " : ", param , " ,AUC = " , auroc, " , prc= ", auprc)
    f.close()      
    
    'Find best in folds and get the mean of param'
    fold_results = np.asanyarray(fold_results)
    mean_param = 0
    for fold_no in xrange(nfold):
        foldi_results = fold_results[fold_results[:,0] == fold_no]
        max_param = foldi_results[np.argmax(foldi_results, axis = 0)[2]]
        mean_param += max_param[1]
        
#         print(max_param)
    mean_param /= nfold
    mean_param_file = outdir + "/" + argsDict["model_type"] + "_mean_param.txt"
    f = open(mean_param_file, 'wt')
    f.write(str(mean_param))
    f.close()
#     print(mean_param)
    
    '============================================='
    '      Plot the cross-validation AUCs'
    '============================================='
    #plot_file = outdir + "/" + argsDict["model_type"] + "_auc_vs_params_crossval.png"
    #plt.figure(figsize=(8, 6), dpi=300)
    #df = pd.DataFrame(fold_results, columns=["fold", "param", "auc", "prc"])
    #fig, ax = plt.subplots()
    #colors = {0:'red', 1:'blue', 2:'green', 3:'black', 4:'orange'}
    #grouped = df.groupby('fold')
    #for key, group in grouped:
    #    group.plot(ax=ax, kind='line', x='param', y='auc', label=key, color=colors[key], xlim = [0, 0.01], ylim = [0.4, 1])
    #plt.savefig(plot_file, dpi=300)
    
    return mean_param

'============================================='
'           Train and Test Method'
'============================================='
def linear_model_simple(train_features, train_labels, test_features, test_labels, param, save_model = False):
    
    logreg = linear_model.LogisticRegression(C=param, penalty='l2')
    logreg.fit(train_features, train_labels)
    y_pred_prob = logreg.predict_proba(test_features)
    
    if save_model:
        coef_outfile = outdir + "/" + argsDict["model_type"] + "_coef_table.txt"
        model_outfile = outdir + "/" + argsDict["model_type"] + "_model.sav"
        train_data_file = outdir + "/" + argsDict["model_type"] + "_x_train"
        coef_vals = logreg.coef_.ravel()
        print(coef_vals.shape , features_header.shape)
        coef_records = np.rec.fromarrays((features_header,coef_vals), names= ('feature', 'coef'))
        coef_df = pd.DataFrame.from_records(coef_records)
        coef_df['abs_coef'] = abs(coef_df['coef'])
        sorted_coef_df = coef_df.sort_values(by = 'abs_coef', ascending = False)
        sorted_coef_df.to_csv(coef_outfile, sep = "\t", columns = ('feature', 'coef'), index = False, header = False)
        pickle.dump(logreg, open(model_outfile, 'wb'))
        np.save(train_data_file, train_features)

    fpr, tpr, thresholds = metrics.roc_curve(test_labels, y_pred_prob[:,1], pos_label=1)
    auroc = metrics.auc(fpr, tpr)
    average_precision = average_precision_score(test_labels, y_pred_prob[:,1])

    return auroc, average_precision

def load_model(model_file):
    # load the model from disk
    loaded_model = pickle.load(open(model_file, 'rb'))

def usage():
    print('train_test_crossval.py -i <input_features.csv> -o <outdir> -n <nfolds> -t <model_type>')

'============================================='
'    Read Inputs'
'============================================='
def readInputs(argv):
    argsDict = {}
    
    try:
        # Reads Array of tuples such as : [("-i", "in.txt"), ("-o", "out.txt")]
        opts, args = getopt.getopt(argv, "hi:o:n:t:", ["infile=", "outfile=", "nfold=", "modelType="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    if len(argv) < 4:
        usage()
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            usage()  
            sys.exit()
        elif opt in ("-i", "--infile"):
            inputFile = arg
        elif opt in ("-o", "--outdir"):
            outputFile = arg
        elif opt in ("-n", "--nfold"):
            nfold = arg
        elif opt in ("-t", "--modelType"):
            model_type = arg
            
    argsDict["infile"] = inputFile
    argsDict["outdir"] = outputFile
    argsDict["nfold"] = nfold
    argsDict["model_type"] = model_type
    return argsDict

'================================================'
'    Main  '
'================================================'
if __name__ == "__main__":
    first_line = True
    argsDict = readInputs(sys.argv[1:])
    with open(argsDict["infile"], 'rt') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            if (first_line):
                features_header = np.asarray(row[0:-1])
                first_line = False
                continue
            all_features.append(np.asarray(row))
        all_features = np.asanyarray(all_features)

    seeds = np.random.randint(120000, size = 3)
    #seeds = [10657, 12421, 541]
    #seed = 3467
    for seed in iter(seeds):
        outdir = argsDict["outdir"] + "/" + str(seed)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
        # split data into 80% train 20% heldout test
        split_data_train_test(seed)
        
        # cross validation and getting mean of best parameters
        mean_param = split_train_cross_val(5)
	#mean_param = 0.000266
        
        # normalize whole train set
    	scaler = preprocessing.StandardScaler(with_mean=True, with_std=True).fit(train_all_features)
        #scaler = preprocessing.MinMaxScaler().fit(train_all_features)
        scaled_train_all = scaler.transform(train_all_features)
        scaled_test_all = scaler.transform(test_all_features)
        auroc, auprc = linear_model_simple(scaled_train_all, train_all_labels, scaled_test_all, test_all_labels, mean_param, save_model = True)
        
        test_outfile = outdir + "/" + argsDict["model_type"] + "_heldout_test.txt"
        f = open(test_outfile, "wt")
        f.write("auroc\tauprc\n")
        f.write(str(auroc) + "\t" + str(auprc) )
        f.close();
        print("heldout test : ", auroc, " prc = ", auprc)
	unscaled_train_file = outdir + "/unscaled_train.npy"
	np.save(unscaled_train_file, train_all_features)
    
    
