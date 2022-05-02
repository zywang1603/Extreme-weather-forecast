path_data = '/Users/mozykhau/University/Extreme event project/Data/'
# Change the above on KNMI server to: /nobackup_1/users/schreurs/project_GAN/

path_project  = '~/Users/mozykhau/PycharmProjects/'
# Change the above on the knmi server to: ~/KNMI_Internship_GANs/

# Global variables that point to the correct directory
dir_rtcor = path_data + '5 min rainfall/RAD_NL25_RAP_5min/' # dataset_radar contains 2019
dir_aart = path_data + 'RAD_NL25_RAC_MFBS_EM_5min/'
# dir_rtcor_npy = path_data + 'dataset_radar_np/'
# dir_aart_npy = path_data + 'dataset_aart_np/'

dir_prep = 'preprocessed/'
dir_rtcor_prep = path_data + dir_prep + 'rad/'
dir_aart_prep = path_data + dir_prep + 'aart/'

dir_labels = path_data + 'rad_rain_labels/'
dir_labels_heavy = path_data + 'rad_heavy_rain_labels/'

prefix_rtcor = 'RAD_NL25_RAP_5min_'
prefix_aart = 'RAD_NL25_RAC_MFBS_EM_5min_'
