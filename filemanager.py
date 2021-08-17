import os, shutil, sys

def GetFileList(dir):
    os.chdir(dir)
    return os.listdir()

def Organize():
    for file in [file for file in os.listdir() if file.find('.') != -1]:
        if not os.path.isdir(file[:file.find('.')]):
            os.mkdir(file[:file.find('.')])
        currDir = os.getcwd()
        os.chdir(os.path.join(os.getcwd(), file[:file.find('.')]))
        shutil.move(os.path.join(currDir, file), os.getcwd())
        print(os.getcwd())
        os.chdir("..")
        print(os.getcwd())


