import subprocess
import time
import datetime
my_file = open("/oak/stanford/groups/leanew1/users/apines/data/gp/DL_list2.txt", "r")
content = my_file.read()
content_list = content. split("\n")
# remove last line (blank)
content_list.pop()
# feed em' in as subjects
subjects = content_list
for x in range(599,len(subjects)):
  newsub = subjects.pop()
  subprocess.call(["sbatch","/oak/stanford/groups/leanew1/users/apines/scripts/gp/FunctionalImages/sbatchABCDProc.sh",newsub])
  time.sleep(700) #wait 20 minutes


# while there are more than 0 subjects left to run
#while len(subjects)>0:
  # grab qsub info, get number of jobs being run
  # qstat = subprocess.check_output(['qstat'],shell=True).decode().split('/bin/python')[0]
  # que = len(qstat.split('\n'))-3
  # if we are using less than 5 job slots (one is occupied by this script)
  # if que < 6:
    # see if it is the weekend, 0, 1, 2, 3, and 4 are weekday, 5 and 6 are weekend
    # weekno = datetime.datetime.today().weekday()
    # see if it is before 9 or after 5 
    # Hour = time.localtime().tm_hour 
    # if weekend OR after 6 PM OR before 8 AM
    # if weekno > 4 or Hour < 8 or Hour > 17 :
    #  newsub = subjects.pop()
      # submit job (if above conditions are met)
    #  subprocess.run(["qsub","-l","h_vmem=11G,s_vmem=10G","qsubMatlab.sh",newsub])
    #  time.sleep(60) #wait a minute

