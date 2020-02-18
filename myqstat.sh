


echo 'job-ID        prior   name       user         state submit/start at     queue                          slots ja-task-ID' 
echo '-----------------------------------------------------------------------------------------------------------------'
qstat -u mtinti -xml | tr '\n' ' ' | sed 's#<job_list[^>]*>#\n#g' | sed 's#<[^>]*>##g' | grep " " | column -t


echo '----'
echo 'tot:' 
qstat -u $1 -xml | grep JB_name  | wc -l

