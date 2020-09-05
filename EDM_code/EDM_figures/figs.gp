#*********************




#*********************
reset
set term pos solid enhanced color font "Times-New-Roman,20" linewidth 2
set out "bL32s100000_stripes.eps"


set xlabel "m_f"
set ylabel "m_p"

unset colorbox

set yrange [-1:1]
set xrange [-.4:.4]

set palette model RGB defined ( -1 'blue', 1 'red' )
p 'bL32s100000/magnetizations.dat' u 2:($3>-1 && $3<1 ? ($3) : 1/0):1 w p pt 7 palette not



#*********************
reset
set term pos solid enhanced color font "Times-New-Roman,20" linewidth 2
set out "bL32s100000_up.eps"

#interval
n=30 #number of intervals
max=.3 #max value
min=-.3 #min value
width=(max-min)/n #interval width

#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5 border lt -1

set xlabel "m_p"
set ylabel "P(m_f=1,m_p)"

set xrange [-.3:.3]

#count and plot
plot "bL32s100000/magnetizations.dat" u (hist($2,width)):($3==1?1.0:0) smooth freq w boxes lc rgb "green" not



#*********************
reset
set term pos solid enhanced color font "Times-New-Roman,20" linewidth 2
set out "bL32s100000_down.eps"

#interval
n=30 #number of intervals
max=.3 #max value
min=-.3 #min value
width=(max-min)/n #interval width

#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5 border lt -1

set xlabel "m_p"
set ylabel "P(m_f=-1,m_p)"

set xrange [-.3:.3]

#count and plot
plot "bL32s100000/magnetizations.dat" u (hist($2,width)):($3==-1?1.0:0) smooth freq w boxes lc rgb "green" not



