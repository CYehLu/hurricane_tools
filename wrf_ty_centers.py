import os


wrfrun_path = '../WRFV3/run/'

command = 'grep ATCF ' + wrfrun_path + 'rsl.error.0000 > centers.txt'
flag = os.system(command)

if flag == 0:
    print('SUCCESSFUL. "centers.txt" is stored at the current folder.')
else:
    print('FAIL')
