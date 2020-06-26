import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
import numpy as np

from chp.core.util import convert_patient_dict_to_dataframe, process_operator


def prepareData(patient_dict, target,
                         include_genes=False,
                         include_gene_variants=False,
                         include_drugs=False,
                         other_info=[]):
    target_name, op_str, val = target
    op = process_operator(op_str)

    y_df_unformat = convert_patient_dict_to_dataframe(patient_dict, other_info=[target_name])

    X_df = convert_patient_dict_to_dataframe(patient_dict,
                                             include_genes=include_genes,
                                             include_gene_variants=include_gene_variants,
                                             include_drugs=include_drugs,
                                             other_info=other_info,
                                             encode_strings=True)

    X = X_df.astype(np.float32).values

    #-- Construct logical test on y data
    conditional = op(y_df_unformat[target_name], val)
    y = conditional.astype(np.float32).values

    return X, y

def randomForestAnalysis(X, y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2, random_state=42)
    clf = RandomForestClassifier(random_state=0)
    clf.fit(X_train, y_train)
    print('Random Forest Score = {}'.format(clf.score(X_test, y_test)))
    return clf

def kNNAnalysis(X,y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2, random_state=42)
    clf = KNeighborsClassifier()
    clf.fit(X_train, y_train)
    print('kNN Score = {}'.format(clf.score(X_test, y_test)))
    return clf

def SVCAnalysis(X,y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2, random_state=42)
    clf = SVC(gamma='auto')
    clf.fit(X_train, y_train)
    print('SVM Score = {}'.format(clf.score(X_test, y_test)))
    return clf

if __name__ == '__main__':
    patient_dict_file ='/home/public/data/ncats/BabelBKBs/collapsedAll/patient_data.pk' 
    with open(patient_dict_file, 'rb') as f_:
        patient_dict = pickle.load(f_)

    target = ('Survival_Time', '>=', 943)
    X, y = prepareData(patient_dict, target, include_genes=True)
    print(X.shape)
    print(y.shape)
    rfc = randomForestAnalysis(X, y)
    knn = kNNAnalysis(X, y)
    svm = SVCAnalysis(X, y)
