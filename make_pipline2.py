# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 11:52:11 2014

@author: mtinti-x
"""
import os

replace_list = eval(open('vars.txt').read())

if __name__ == '__main__':
    run_all_content = ''
    fast_q_files = []
    for dictionary in replace_list:
        fast_q_files.append(dictionary['base_fastq'])
        sh_script_name = dictionary['base_fastq']+'.sh'
        run_all_content+='chmod +x '+sh_script_name+'\n'
        run_all_content+='qsub '+sh_script_name+'\n'
        run_all_content+='mv '+sh_script_name+' '+dictionary['experiment']+'\n'
        
        template_file = open(os.path.join('templates','_template_mapping_multiple2.sh')).read()
        sh_script_content = template_file.replace(
            '{_genome_}',dictionary['g_version']).replace(
            '{fastq}',dictionary['base_fastq']).replace(
            '{experiment}',dictionary['experiment']).replace(
            '{library}',dictionary['library'])
        open(sh_script_name,'w').write(sh_script_content)
        

    open('run_all_'+dictionary['experiment']+'.sh','w').write(run_all_content)

    
