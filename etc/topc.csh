
set DATE = `date +%H%M`
echo $DATE

ps -e -o user -o time -o comm | sort -k 2 >! $1/topcpu.$DATE


