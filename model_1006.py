# -*- coding: utf-8 -*-
"""
    批量链接数据以及修改数据名 (Python 版改名，不依赖系统 rename 命令)
"""
import os
import subprocess
import glob

def CopyFile(origfile_filedir, simfile_rootfiledir, origfile_filename, simfile_filenames, origfile_year, simfile_year):
    """
        批量链接数据以及修改数据名
    """
    for simfile_filename in simfile_filenames:
        simfile_filedir = os.path.join(simfile_rootfiledir, simfile_filename)
        if not os.path.exists(simfile_filedir):
            os.makedirs(simfile_filedir)
        
        # 1. 链接 landdata
        origfile_filepath = os.path.join(origfile_filedir, "landdata")
        simfile_filepath = os.path.join(simfile_filedir, "landdata")
        
        # 删除旧的 landdata 链接或目录，防止多层嵌套
        if os.path.islink(simfile_filepath) or os.path.exists(simfile_filepath):
            subprocess.run(f"rm -rf {simfile_filepath}", shell=True, check=True)

        # 创建新的软链接
        p = subprocess.Popen("ln -sf {} {}".format(origfile_filepath, simfile_filepath), shell=True)
        p.wait()
        
        # 2. 拷贝 restart
        simfile_restart = os.path.join(simfile_filedir, "restart")
        if not os.path.exists(simfile_restart):
            os.makedirs(simfile_restart)

        # 清理旧的 restart 目录，防止出现多个年份文件夹
        if os.path.exists(simfile_restart):
            subprocess.run(f"rm -rf {simfile_restart}", shell=True, check=True)
        os.makedirs(simfile_restart)
        
        # 年份目录
        orig_restart_dir = os.path.join(origfile_filedir, "restart", "{}-001-00000".format(origfile_year))
        new_year_dir = os.path.join(simfile_restart, "{}-001-00000".format(simfile_year))
        p = subprocess.Popen("cp -r {} {}".format(orig_restart_dir, new_year_dir), shell=True)
        p.wait()

        # const 目录
        orig_const_dir = os.path.join(origfile_filedir, "restart", "const")
        p = subprocess.Popen("cp -r {} {}/".format(orig_const_dir, simfile_restart), shell=True)
        p.wait()

        # 3. Python 改名
        # 3.1 改 const 目录下的文件名
        const_dir = os.path.join(simfile_restart, "const")
        for f in glob.glob(os.path.join(const_dir, "*")):
            new_name = f.replace(origfile_filename, simfile_filename)
            if new_name != f:
                os.rename(f, new_name)

        # 3.2 改年份目录名
        old_year_dir = os.path.join(simfile_restart, "{}-001-00000".format(origfile_year))
        new_year_dir = os.path.join(simfile_restart, "{}-001-00000".format(simfile_year))
        if os.path.exists(old_year_dir):
            os.rename(old_year_dir, new_year_dir)

        # 3.3 改年份目录里的文件名
        if os.path.exists(new_year_dir):
            for f in glob.glob(os.path.join(new_year_dir, "*")):
                new_name = f.replace(origfile_filename, simfile_filename).replace(str(origfile_year), str(simfile_year))
                if new_name != f:
                    os.rename(f, new_name)

        print("CopyFile: {} is done!".format(simfile_filename))


if __name__ == "__main__":
    origfile_filedir = "/stu02/zhangsy24/model/CoLMcases/TEST20251201/US-Ton_mkinidata_spinup_1/"
    simfile_rootfiledir = "/stu02/zhangsy24/model/CoLMcases/TEST20251201/"
    simfile_filenames = ["US-Ton_4_1land8_nospin_T1"]
    origfile_filename = "US-Ton_mkinidata_spinup_1"
    origfile_year = 2020
    simfile_year = 1980

    CopyFile(origfile_filedir, simfile_rootfiledir, origfile_filename, simfile_filenames, origfile_year, simfile_year)
