#! /bin/tcsh -f

set name = `ls $1 |awk -F_ '{print $1}'`
echo $name
touch $name
set pippi = `wc $1 |  awk '{print $1}'`
set countino = 0

while ( $countino < $pippi )

   set line = `head -$countino $1 |  tail -1`
 
   set value = `echo $line | awk -F/ '{print $6}' |  awk -F_ '{print $2}'`
   set BRBR = `echo $line | awk -F= '{print $2}' | awk -F'+-' '{print $1}'`
   set errBRBR = `echo $line | awk -F= '{print $2}' | awk -F'+-' '{print $2}' | awk -F' ' '{print $2}' |  awk -F'(stat)' '{print $1}'`
   set errstatBRBR = `echo $line | awk -F= '{print $2}' | awk -F'+-' '{print $3}' | awk -F' ' '{print $2}' | awk -F'(stat)' '{print $1}'`
    
    echo $value $BRBR $errBRBR $errstatBRBR >> $name

   @ countino ++

end
