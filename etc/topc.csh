
set DATE = `date +%H%M`

ps -ef -o user -o time -o comm | sort -k 2 >! $1/topcpu.$DATE

uptime >> $1/topcpu.$DATE
