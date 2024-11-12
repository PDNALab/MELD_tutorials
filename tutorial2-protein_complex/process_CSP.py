final=[3,6,9,20,23,24,25,33,35,36,37,38,41,42,43,44,45,46,47]

ET=str('SYDEKRQLSLDINRLPGEKLGRVVHIIQSREPSLRDSNPDEIEIDFETLKPTTLRELERYVKSCLQKK') #TP
minprot=str('MKEEEIEEMIEKAKEELRKRYPEAKEVFLSFTYEVNGKLKIKLTRFDPSMSLEEVEERIEEEVKRLLKEADSIEIRVHTTV')
seq=ET+minprot

#all

hairpin1=25
hairpin2=47

file=open('prot_pep_all_fixed.dat','w')
for i in final:
    for j in range(69,150): # add 1
        if j < hairpin1+69 or j > hairpin2+69:
            continue
        if seq[i-1] =='G' and seq[j-1] != 'G':
            file.write( "{} CA {} CB {}\n".format(i, j, 0.8 ))
            file.write("\n")
        if seq[i-1] !='G' and seq[j-1] =='G':
            file.write( "{} CB {} CA {}\n".format(i, j, 0.8 ))
            file.write("\n")
        if seq[i-1] =='G' and seq[j-1] =='G':
            file.write( "{} CA {} CA {}\n".format(i, j, 0.8 ))
            file.write("\n")
        if seq[i-1] !='G' and seq[j-1] !='G':
            file.write( "{} CB {} CB {}\n".format(i, j, 0.8 ))
            file.write("\n")
    
file.close()

