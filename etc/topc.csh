
set DATE = `date +%H%M`

ps -e -o user -o time -o comm | sort -k 2 >! $1/topcpu.$DATE


