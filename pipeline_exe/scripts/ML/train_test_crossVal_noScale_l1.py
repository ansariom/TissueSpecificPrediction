import sys, getopt
import numpy as np
import csv
from sklearn import preprocessing
from numpy import float32, int32
from sklearn.preprocessing.data import StandardScaler
from sklearn import metrics
from sklearn import linear_model
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import average_precision_score
from os.path import os
import pickle
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split

max_param = 3
all_features = []
test_all_features = []
train_all_labels = []
test_all_labels = []
train_all_features = []
features_header = []
batch_index = 0
learning_rate = 0.1
training_epochs = 200
batch_size = 100
display_step = 50
weight_decays = np.arange(0.01, max_param ,0.1) 
fold = 0
fold_results = []
outdir = "-"

'Split data into train and helout test'
def split_data_train_test(seed):
    global train_all_features, test_all_features, train_all_labels, test_all_labels
    print(seed)

    #np.random.seed(seed)
    all_features_ = np.asarray(all_features[:, 1:-1], float32)
    all_labels_ = np.asarray(all_features[:, -1],int32).reshape(len(all_features), 1)
    
    np.random.seed(seed)
    train_all_features, test_all_features, train_all_labels, test_all_labels = train_test_split(all_features_, all_labels_, test_size=0.2)

    train_tss_list = outdir + "/" + argsDict["model_type"] + "_train_tss_ids.txt"
    test_tss_list = outdir + "/" + argsDict["model_type"] + "_test_tss_ids.txt"
    np.savetxt(train_tss_list, train_all_features[:,0], delimiter = "\t", fmt="%s")
    np.savetxt(test_tss_list, test_all_features[:,0], delimiter = "\t" , fmt = "%s")
    
def split_train_cross_val(nfold):
    'split train into a pie with 5 section, 4 train 1 validation'
    param_file = outdir + "/" + argsDict["model_type"] + "_crossval_params_performace.txt"
    f = open(param_file, "wt")
   
    fold_results = []
    fold_tprs = []
    fold_fprs = []
    fold_aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    
    fold_num = 0
    skf = StratifiedKFold(n_splits=5)
    skf.get_n_splits(train_all_features, train_all_labels)
    
    i = 0
    for train_index, test_index in skf.split(train_all_features, train_all_labels):
	train_fold_features, test_fold_features = train_all_features[train_index], train_all_features[test_index]
        train_fold_labels, test_fold_labels = train_all_labels[train_index], train_all_labels[test_index]
        
	tprs = []
        fprs = []
    	aucs = []
	for param in weight_decays:
	    auroc, auprc, tpr, fpr = linear_model_simple(train_fold_features, train_fold_labels, test_fold_features, test_fold_labels, param)
	    fold_results.append(np.array([fold_num, param, auroc, auprc]))     
	    f.write(str(fold_num) + "\t" + str(param) + "\t" + str(auroc) + "\t" + str(auprc) + "\n")
	    
	    #fprs.append(fpr)
	    tprs.append(np.interp(mean_fpr, fpr, tpr))
	    tprs[-1][0] = 0.0
	    aucs.append(auroc)
            # plt.plot(fpr, tpr, lw=1, alpha=0.3, label='ROC fold %d (AUC = %0.2f)' % (i, auroc))
    	print(fold_results)
        fold_tprs.append(np.mean(tprs, axis=0))
	fold_fprs.append(mean_fpr)
	#fold_fprs.append(np.mean(fprs, axis=0))
	fold_aucs.append(np.max(aucs, axis=0))
	plt.plot(mean_fpr, np.mean(tprs, axis=0), lw=1, alpha=0.3, label='ROC fold %d (AUC = %0.2f)' % (fold_num, np.max(aucs, axis=0)))
    	fold_num += 1

    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)
    mean_tpr = np.mean(fold_tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = metrics.auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    plt.plot(mean_fpr, mean_tpr, color='b',
         label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
         lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                 label=r'$\pm$ 1 std. dev.')

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic (ROC) for 5-fold CV')
    plt.legend(loc="lower right")
    plot_file = outdir + "/" + argsDict["model_type"] + "_auc_vs_params_crossval.png"
    #plt.figure(figsize=(8, 6), dpi=300)
    #fig, ax = plt.subplots()
    plt.savefig(plot_file, dpi=300)
    plt.close()

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
    #    group.plot(ax=ax, kind='line', x='param', y='auc', label=key, color=colors[key], xlim = [0, 0.1], ylim = [0.4, 1])
    #plt.savefig(plot_file, dpi=300)
    
    return mean_param

'============================================='
'           Train and Test Method'
'============================================='
def linear_model_simple(train_features, train_labels, test_features, test_labels, param, save_model = False):
    logreg = linear_model.LogisticRegression(C=param, penalty='l1')
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
    return auroc, average_precision, tpr, fpr

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

def read_input(argsDict):
    first_line = True
    all_data = []
    with open(argsDict["infile"], 'rt') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            if (first_line):
                features_header = np.asarray(row[0:-1])
                first_line = False
                continue
            all_data.append(np.asarray(row))
    all_features = np.asanyarray(all_data)
    return all_features, features_header
   

'================================================'
'    Main  '
'================================================'
if __name__ == "__main__":
    argsDict = readInputs(sys.argv[1:])
    all_features, features_header = read_input(argsDict)
    seeds = np.random.randint(1200000, size = 30)
    #seeds = [10657, 12421, 541]
    #seeds = [107244, 30222911, 593221000, 67, 230] 
#    seeds = [107244] 
    for seed in iter(seeds):
        outdir = argsDict["outdir"] + "/" + str(seed)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
        # split data into 80% train 20% heldout test
        split_data_train_test(seed)
        
        # cross validation and getting mean of best parameters
        mean_param = split_train_cross_val(5)
	#mean_param = 0.08001
        
        # normalize whole train set
#    	scaler = preprocessing.StandardScaler(with_mean=True, with_std=True).fit(train_all_features)
#        scaler = preprocessing.MinMaxScaler().fit(train_all_features)
#        scaled_train_all = scaler.transform(train_all_features)
#        scaled_test_all = scaler.transform(test_all_features)
	scaled_train_all = train_all_features
	scaled_test_all = test_all_features
        auroc, auprc, tpr, fpr = linear_model_simple(scaled_train_all, train_all_labels, scaled_test_all, test_all_labels, mean_param, save_model = True)
        
        test_outfile = outdir + "/" + argsDict["model_type"] + "_heldout_test.txt"
        f = open(test_outfile, "wt")
        f.write("auroc\tauprc\n")
        f.write(str(auroc) + "\t" + str(auprc) )
        f.close();
        print("heldout test : ", auroc, " prc = ", auprc)
	unscaled_train_file = outdir + "/unscaled_train.npy"
	np.save(unscaled_train_file, train_all_features)
    
    
