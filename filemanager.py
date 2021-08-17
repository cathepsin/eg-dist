import os, shutil, sys, ntpath

def GetFileList(dir):
    os.chdir(dir)
    return os.listdir()

def Organize():
    for file in [file for file in os.listdir() if file.rfind('.') != -1]:
        if not os.path.isdir(file[:file.rfind('.')]):
            os.mkdir(file[:file.rfind('.')])
        currDir = os.getcwd()
        os.chdir(os.path.join(os.getcwd(), file[:file.rfind('.')]))
        shutil.move(os.path.join(currDir, file), os.getcwd())
        print(os.getcwd())
        os.chdir("..")
        print(os.getcwd())

def TwoFiles(f1, f2):
    if not os.path.isdir(f1[:f1.rfind(".")] + " eg-dist"):
        os.mkdir(f1[:f1.rfind(".")] + " eg-dist")
    os.chdir(os.path.join(os.getcwd(), f1[:f1.rfind(".")] + " eg-dist"))
    print(os.path.basename(f1), os.path.basename(f2))
    if not os.path.isfile(os.path.basename(f1)):
        shutil.copy(f1, os.getcwd())
    if not os.path.isfile(os.path.basename(f2)):
        shutil.copy(f2, os.getcwd())
    os.chdir("..")
    return f1[:f1.rfind(".")] + " eg-dist"


