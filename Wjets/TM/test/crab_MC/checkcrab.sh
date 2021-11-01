for i in `ls crab_* -d`
do
echo $i
#crab status -d $i
sleep 5
crab kill -d $i
done
