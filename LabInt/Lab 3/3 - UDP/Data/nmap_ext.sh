cat packets.txt | tr '→' ':' | tr -s ':'| cut -d ':' -f 2 | cut -d ' ' -f 2 > ports.txt
cat packets.txt | tr -s ' ' | cut -d ' ' -f 3 > time.txt
paste time.txt ports.txt > Data.txt
rm time.txt
rm ports.txt
