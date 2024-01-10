# _*_ coding: utf-8 _*_
# author: whitebird
# date: 2020/9/2
#-----------------------------------------------------------------------------#

import argparse, os, re
import pandas as pd
from function_basic import *


def runCellranger(project_name, project_species):
    # project_species = args.species
    config_dic = readconfig()
    cellranger_data_path = config_dic['cellranger_data']
    cellranger_result_path = config_dic['cellranger_result']
    cellranger_data_ftp = config_dic['cellranger_data_ftp']
    cellranger_ftp = config_dic['cellranger_ftp']
    project_path = config_dic['project_path']
    cellranger_lab = config_dic['cellranger_lab']

    project_cellranger_data_dir = os.path.join(cellranger_data_path, project_name)
    project_cellranger_result_dir = os.path.join(cellranger_result_path, project_name)
    project_analysis_path = os.path.join(project_path, project_name)
    project_cellranger_data_ftp = os.path.join(cellranger_data_ftp, project_name)
    project_cellranger_ftp = os.path.join(cellranger_ftp, project_name)
    project_cellranger_lab = os.path.join(cellranger_lab, project_name)

    for sample in os.listdir(project_cellranger_data_dir):
        if os.path.isdir(os.path.join(project_cellranger_data_dir, sample)):
            project_sample_dir = os.path.join(project_cellranger_data_dir, sample)
            project_sample_file = [item for item in os.listdir(project_sample_dir) if item.endswith('.gz')]

            # 代码进行调整，物种信息为：human，mouse，human_mouse；如果为human_mouse，跑human,mouse,human_mouse三种；
            project_species_list = []
            if project_species.lower() == "human":
                project_species_list = ["human"]
            elif project_species.lower() == "mouse":
                project_species_list = ["mouse"]
            elif project_species.lower() == "hcmv":
                project_species_list = ["hcmv"]
            elif project_species.lower() == "human_hcmv":
                project_species_list = ["human_hcmv"]    
            elif project_species.lower() == "human_mouse":
                project_species_list = ["human", "mouse", "human_mouse"]
            # sample_name转换成sample_id
            if df1.loc[(df1["pro_code"] == project_name.strip()) & (df1["sample_id"] == sample.strip())].sample_name.count()==0:
                continue
            # sample_name = df1.loc[(df1["pro_code"] == project_name.strip()) & (df1["sample_id"] == sample.strip())].sample_name.max()
            # for item in project_sample_file:
            #     os.rename(os.path.join(project_cellranger_data_dir, sample, item),
            #               os.path.join(project_cellranger_data_dir, sample, item.replace(sample_name, sample)))
                
            seq_type = df1.loc[
                (df1["pro_code"] == project_name.strip()) & (df1["sample_id"] == sample.strip())].seq_type.max()
            # rawdata备份目录
            sample_cellranger_data_ftp = os.path.join(project_cellranger_data_ftp, 'Rawdata', seq_type, sample)
            # 样本数据的ftp备份
            if not os.path.exists(sample_cellranger_data_ftp):
                os.makedirs(sample_cellranger_data_ftp, 0o775)

            # part02: 生成bash脚本，投递任务
            seq_type = df1.loc[
                (df1["pro_code"] == project_name.strip()) & (df1["sample_id"] == sample.strip())].seq_type.max()
            print(sample, project_name, seq_type)
            for project_species_i in project_species_list:
                if len(project_species_list) == 1:
                    dir_type = seq_type.lower()
                else:
                    dir_type = project_species_i.lower()
                if os.path.exists(os.path.join(project_cellranger_result_dir, dir_type, sample, 'outs')):
                    continue
                sample_sh = project_analysis_path + f"/run_CellRanger_{sample}_{project_species_i}.sh"
                print(sample, project_name, project_species_i)
                with open(sample_sh, "w") as sh_handel:
                    sh_handel.write("#!/bin/bash\n")
                    sh_handel.write("#$ -S /bin/bash\n")
                    sh_handel.write("#$ -V\n")
                    sh_handel.write("#$ -cwd\n")
                    sh_handel.write("#$ -l vf=50G\n")
                    # sh_handel.write(f"#$-N {sample}_CellRanger\n")
                    # sh_handel.write("#$-q all.q@fat01,all.q@cu01,all.q@cu02,all.q@cu03,all.q@cu04,all.q@cu05,all.q@cu06,all.q@cu07,all.q@cu08\n")
                    # sh_handel.write("#$-l vf=36g,p=4\n")
                    sh_handel.write(f"cellranger='{config_dic['cellranger']}'\n")
                    if seq_type in ['tcr', 'bcr', 'TCR', 'BCR']:
                        if project_species_i.lower() == "human":
                            sh_handel.write(f"transcriptome='{config_dic['vdj_transcriptome_human']}'\n")
                        elif project_species_i.lower() == "mouse":
                            sh_handel.write(f"transcriptome='{config_dic['vdj_transcriptome_mouse']}'\n")
                        elif project_species_i.lower() == "human_mouse":
                            sh_handel.write(f"transcriptome='{config_dic['vdj_transcriptome_human_mouse']}'\n")
                    else:
                        if project_species_i.lower() == "human":
                            sh_handel.write(f"transcriptome='{config_dic['transcriptome_human']}'\n")
                        elif project_species_i.lower() == "mouse":
                            sh_handel.write(f"transcriptome='{config_dic['transcriptome_mouse']}'\n")
                        elif project_species_i.lower() == "human_mouse":
                            sh_handel.write(f"transcriptome='{config_dic['transcriptome_human_mouse']}'\n")
                        elif project_species_i.lower() == "aaegl5":
                            sh_handel.write(f"transcriptome='{config_dic['transcriptome_AaegL5']}'\n")
                        elif project_species_i.lower() == "hcmv":
                            sh_handel.write(f"transcriptome='{config_dic['transcriptome_HCMV']}'\n")
                        elif project_species_i.lower() == "human_hcmv":
                            sh_handel.write(f"transcriptome='{config_dic['transcriptome_human_HCMV']}'\n")       
                    sh_handel.write(f"cd {project_sample_dir};\n")
                    sh_handel.write(
                        f"if [ ! -f '/{sample}_md5sum.txt' ];then md5sum ./*fastq.gz >{sample}_md5sum.txt; fi\n")
                    sh_handel.write(f'cp *.gz  {sample_cellranger_data_ftp};\n')
                    sh_handel.write(f"mkdir -p {project_cellranger_result_dir}/{dir_type};\n")
                    sh_handel.write(f"cd {project_cellranger_result_dir}/{dir_type};\n")
                    sh_handel.write(f"rm -rf {sample};\n")
                    if seq_type in ['tcr', 'TCR']:
                        sh_handel.write(
                            f"$cellranger vdj --id={sample} --fastqs={project_sample_dir} --reference=$transcriptome --sample={sample} --chain  TR;\n")
                    elif seq_type in ['bcr', 'BCR']:
                        sh_handel.write(
                            f"$cellranger vdj --id={sample} --fastqs={project_sample_dir} --reference=$transcriptome --sample={sample} --chain  IG;\n")
                    elif seq_type in ['rna_seq']:
                        sh_handel.write(
                            f"$cellranger count --id={sample} --fastqs={project_sample_dir} --transcriptome=$transcriptome --sample={sample} --localcores 12 --localmem 128;\n")
                    elif seq_type in ['snrna_seq']:
                        sh_handel.write(
                            f"$cellranger count --id={sample} --fastqs={project_sample_dir} --transcriptome=$transcriptome --sample={sample} --localcores 8 --localmem 64 --include-introns;\n")
                    else:
                        sh_handel.write(
                            f"$cellranger count --id={sample} --fastqs={project_sample_dir} --transcriptome=$transcriptome --sample={sample} --localcores 8 --localmem 64 --include-introns --chemistry=ARC-v1;\n")

                    # 客户版本共享文件 (02_cellranger_result)
                    sh_handel.write(f"mkdir -p {project_cellranger_ftp}/{dir_type}/{sample}/;\n")
                    # 内部版本共享文件 (02_cellranger_lab)
                    sh_handel.write(f"mkdir -p {project_cellranger_lab}/{dir_type};\n")
                    # 内部版本共享文件拷贝网页文件
                    sh_handel.write(
                        f"cp ./{sample}/outs/web_summary.html {project_cellranger_lab}/{dir_type}/{sample}.html;\n")
                    # 02_cellranger_result
                    if seq_type in ['tcr', 'bcr', 'TCR', 'BCR']:
                        sh_handel.write(
                            f"cp ./{sample}/outs/web_summary.html {project_cellranger_ftp}/{dir_type}/{sample};\n")
                        sh_handel.write(f"cp ./{sample}/outs/vloupe.vloupe {project_cellranger_ftp}/{dir_type}/{sample};\n")
                        # sh_handel.write(f"python /home/project/Single_cells_2020/scrna_project/pipeline/py_script/edit_html.py -frompath {project_cellranger_internal}/{seq_type}/{sample} -topath {project_cellranger_ftp}/{seq_type}/{sample};\n")
                        sh_handel.write(f"cp ./{sample}/outs/*.csv {project_cellranger_ftp}/{dir_type}/{sample};\n")
                    else:
                        sh_handel.write(
                            f"cp -r ./{sample}/outs/filtered_feature_bc_matrix {project_cellranger_ftp}/{dir_type}/{sample};\n")
                        sh_handel.write(
                            f"cp ./{sample}/outs/metrics_summary.csv {project_cellranger_ftp}/{dir_type}/{sample};\n")
                        sh_handel.write(
                            f"cp ./{sample}/outs/web_summary.html {project_cellranger_ftp}/{dir_type}/{sample};\n")
                        # sh_handel.write(f"python /home/project/Single_cells_2020/scrna_project/pipeline/py_script/edit_html.py -frompath {project_cellranger_internal}/{seq_type}/{sample} -topath {project_cellranger_ftp}/{seq_type}/{sample};\n")
                    # 项目文件 (scrna_project)
                    sh_handel.write(f"mkdir -p {project_analysis_path}/data/{dir_type}/{sample};\n")
                    if seq_type in ['tcr', 'bcr', 'TCR', 'BCR']:
                        sh_handel.write(
                            f"cp ./{sample}/outs/web_summary.html {project_analysis_path}/data/{dir_type}/{sample};\n")
                        sh_handel.write(f"cp ./{sample}/outs/*.csv {project_analysis_path}/data/{dir_type}/{sample};\n")
                        sh_handel.write(
                            f"cp ./{sample}/outs/vloupe.vloupe {project_analysis_path}/data/{dir_type}/{sample};\n")
                    else:
                        sh_handel.write(
                            f"cp -r ./{sample}/outs/filtered_feature_bc_matrix {project_analysis_path}/data/{dir_type}/{sample};\n")
                        sh_handel.write(
                            f"cp ./{sample}/outs/metrics_summary.csv {project_analysis_path}/data/{dir_type}/{sample};\n")
                        sh_handel.write(
                            f"cp ./{sample}/outs/web_summary.html {project_analysis_path}/data/{dir_type}/{sample};\n")
                    # 为项目的rawdata生成MD5值
                    sh_handel.write(f"cd {project_cellranger_data_ftp};\n")
                    sh_handel.write(f"md5sum ./*/*/*/*fastq.gz >MD5.txt;\n")
                # 项目目录作为当前目录
                os.chdir(project_analysis_path)
                # os.system(f"qsub run_CellRanger.sh")
                # os.system(f"qsub -N Ce_{sample} -q all.q@fat01,all.q@cu01,all.q@cu03,all.q@cu04,all.q@cu05,all.q@cu06,all.q@cu07 run_CellRanger_{sample}.sh")
                os.system(f"qsub -N Ce_{sample} -q batch@fat01,batch@cu01,batch@cu03,batch@cu04,batch@cu05,batch@cu06 run_CellRanger_{sample}_{project_species_i}.sh")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-project', '--project', type=str, required=True)
    #parser.add_argument('-species', '--species', type=str, required=True, default="human")

    args = parser.parse_args()
    project_name = args.project.strip()
    config_dic = readconfig()
    filename = config_dic['project_filename']
    df1 = pd.DataFrame(pd.read_excel(filename, sheet_name='测序量统计'))
    df1.columns = ['pro_name', 'pro_unit', 'pro_customer', 'pro_code', 'sample_name', 'sample_id', 'seq_type',
                   'species', 'organization', 'is_done']
    df1 = df1[(df1[u'is_done'] != 'T')]
    print(df1.head())

    df1["pro_code"] = df1["pro_code"].astype('str')
    df1["sample_name"] = df1["sample_name"].astype('str')
    df1["sample_id"] = df1["sample_id"].astype('str')
    project_species = df1.loc[df1["pro_code"] == project_name.strip()].species.max()
    print('project_name:', project_name, 'project_species:', project_species)

    runCellranger(project_name, project_species)

