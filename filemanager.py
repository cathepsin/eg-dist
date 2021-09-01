import os, shutil

#Return list of files in directory. Also navigates to directory
def GetFileList(dir):
    os.chdir(dir)
    return os.listdir()

#Organize all given files by name into new directories
def Organize():
    for file in [file for file in os.listdir() if file.rfind('.') != -1]:
        if not os.path.isdir(file[:file.find('.')]):
            os.mkdir(file[:file.find('.')])
        currDir = os.getcwd()
        os.chdir(os.path.join(os.getcwd(), file[:file.find('.')]))
        shutil.move(os.path.join(currDir, file), os.getcwd())
        os.chdir("..")

#If provided two files using the -file flag
def TwoFiles(f1, f2):
    if not os.path.isdir(f1[:f1.rfind(".")] + " eg-dist"):
        os.mkdir(f1[:f1.rfind(".")] + " eg-dist")
    os.chdir(os.path.join(os.getcwd(), f1[:f1.rfind(".")] + " eg-dist"))
    if not os.path.isfile(os.path.basename(f1)):
        shutil.copy(f1, os.getcwd())
    if not os.path.isfile(os.path.basename(f2)):
        shutil.copy(f2, os.getcwd())
    os.chdir("..")
    return f1[:f1.rfind(".")] + " eg-dist"