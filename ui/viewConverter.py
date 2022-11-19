from os import listdir, getcwd, system
path = getcwd()
for file in (file for file in listdir() if file[-3:] == '.ui'):
    system(f'pyuic5 '+path+'\\'+file+' -o '+path+'\\'+file[:-3]+'.py')