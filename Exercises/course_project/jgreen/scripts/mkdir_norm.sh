declare -a StringArray3=("CA_J06" "CA_J08" "CA_J11" "CA_J18G" "CA_J18" "CASE_J03" "CASE_J09" "CASE_J12" "CASE_J13G" "CASE_J13" "CON_J02" "CON_J05" "CON_J10" "CON_J15" "IS_01" "IS_02" "IS_03_" "IS_04" "SE_J01" "SE_J04" "SE_J07" "SE_J14G" "SE_J14")
for i in "${StringArray3[@]}"
do
mkdir $i
mv ${i}*.fq ${i}/
done

for i in "${StringArray3[@]}"
do