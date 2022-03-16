'''
#This is the original config.py
path_data = '/ceph/knmimo/'
# Change the above on KNMI server to: /nobackup_1/users/schreurs/project_GAN/

path_project  = '~/Documents/KNMI_Internship_GANs/'
# Change the above on the knmi server to: ~/KNMI_Internship_GANs/

# Global variables that point to the correct directory
dir_rtcor = path_data + 'dataset_rtcor/' # dataset_radar contains 2019
dir_aart = path_data + 'dataset_aart/'

dir_rtcor_npy = path_data + 'dataset_radar_np/'
dir_aart_npy = path_data + 'dataset_aart_np/'

dir_prep = 'preprocessed/'
dir_rtcor_prep = path_data + dir_prep + 'rtcor/'
dir_aart_prep = path_data + dir_prep + 'aart/'

dir_labels = path_data + 'rtcor_rain_labels/'
dir_labels_heavy = path_data + 'rtcor_heavy_rain_labels/'

prefix_rtcor = 'RAD_NL25_RAC_5M_'
prefix_aart = 'RAD_NL25_RAC_MFBS_EM_5min_'

'''

#The following part is written by Zhiyi Wang on March 14th 2022
path_data = 'E:/thesis_datasets/'
# Change the above on KNMI server to: /nobackup_1/users/schreurs/project_GAN/

path_project  = 'E:/thesis_code/Large_Sample_Nowcasting_Evaluation/'
# Change the above on the knmi server to: ~/KNMI_Internship_GANs/

# Global variables that point to the correct directory
dir_rtcor = path_data  # dataset_radar contains 2019
dir_aart = path_data + 'RADNL_CLIM_EM_MFBSNL25_05m_20081231T235500_20091231T235500_0002/RAD_NL25_RAC_MFBS_EM_5min/'

dir_rtcor_npy = path_data + 'dataset_radar_np/'
dir_aart_npy = path_data + 'dataset_aart_np/'


dir_prep = 'preprocessed/'
dir_rtcor_prep = path_data + dir_prep + 'rtcor/'
dir_aart_prep = path_data + dir_prep + 'aart/'

dir_labels = path_data + 'rtcor_rain_labels/'
#dir_labels_heavy = path_data + 'rtcor_heavy_rain_labels/'
dir_labels_heavy = dir_labels

prefix_rtcor = 'RAD_NL25_RAP_5min_'
prefix_aart = 'RAD_NL25_RAC_MFBS_EM_5min_'
