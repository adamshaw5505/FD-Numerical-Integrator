FILES = system("ls -1 *.dat") 
LABEL = system("ls -1 *dat")
set logscale y
plot for [i=1:lines(FILES)] word(FILES,i) u 1:2 title word(LABEL,i) w l 
pause -1
