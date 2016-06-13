for i in out/*out; do 
  printf "%s " $i
  head -4 $i | tail -n 1 
done 
