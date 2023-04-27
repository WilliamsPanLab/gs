# set list to copy to new list
cp /oak/stanford/groups/leanew1/users/apines/data/gp/DL_list.txt /oak/stanford/groups/leanew1/users/apines/data/gp/DL_list2.txt

for subject in $(cat /oak/stanford/groups/leanew1/users/apines/data/gp/DL_list2.txt); do
    if [ -f "/oak/stanford/groups/leanew1/users/apines/data/gp/PropFeats/$subject/Prop_Feats.csv" ]; then
        sed -i "/$subject/d" /oak/stanford/groups/leanew1/users/apines/data/gp/DL_list2.txt
    fi
done

