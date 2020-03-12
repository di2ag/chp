import matplotlib.pyplot as plt
import pandas as pd
import os

DIR = '/home/public/data/ncats/ValidationTests'

test_names = {'Tamoxifen': {300: 'crossValid-results-tamoxifen-300.csv',
                           1000: 'crossValid-results-tamoxifen-1000.csv'},
              'HERCEPTIN': {300: 'crossValid-results-HERCEPTIN.csv',
                           1000: 'crossValid-results-HERCEPTINSurv1000.csv'},
              'CYCLOPHOSPHAMIDE': {300: 'crossValid-results-CYCLOPHOSPHAMIDE.csv',
                           1000: 'crossValid-results-CYCLOPHOSPHAMIDESurv1000.csv'},}

test_names_1 = {'No Drugs': {300: 'crossValid-results-noDrugs-300.csv',
                           1000: 'crossValid-results-noDrugs-1000.csv'}}
fig, ax = plt.subplots(2)

ax[0].set_title('Prob Survival >= 1YR')
ax[0].set_xlabel('Age at Diagnosis (days)')
ax[0].set_ylabel('Probability')
for drug, result_dict in test_names.items():
    df = pd.read_csv(os.path.join(DIR, result_dict[300]))
    ax[0].plot(df.loc[:,'X'], df.loc[:,"('Updates', 'Survival_Time >= 300', 'True')"], label=drug)
    ax[0].legend()

ax[1].set_title('Prob Survival >= 3YR')
ax[1].set_xlabel('Age at Diagnosis (days)')
ax[1].set_ylabel('Probability')
for drug, result_dict in test_names.items():
    df = pd.read_csv(os.path.join(DIR, result_dict[1000]))
    df.replace(-1,0, inplace=True)
    ax[1].plot(df.loc[:,'X'], df.loc[:,"('Updates', 'Survival_Time >= 1000', 'True')"], label=drug)
    ax[1].legend()
plt.show()


fig, ax = plt.subplots(1)
ax.set_title('Prob Survival with No Drugs.')
ax.set_xlabel('Age at Diagnosis (days)')
ax.set_ylabel('Probability')
for _, result_dict in test_names_1.items():
    for days, file_ in result_dict.items():
        df = pd.read_csv(os.path.join(DIR, file_))
        ax.plot(df.loc[:,'X'], df.loc[:,"('Updates', 'Survival_Time >= {}', 'True')".format(days)], label='Survival >= {}'.format(days))
        ax.legend()

plt.show()
