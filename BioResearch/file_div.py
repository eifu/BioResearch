 

writefile2= open("wig/MM8_mES_WT_shRNAGFP_shRNASetDB1_normalized_H3K9me3_two.WIG","w")
writefile3= open("wig/MM8_mES_WT_shRNAGFP_shRNASetDB1_normalized_H3K9me3_three.WIG","w")

        
with open("wig/MM8_mES_WT_shRNAGFP_shRNASetDB1_normalized_H3K9me3.WIG","r") as f:
    firstline = f.readline()
    secondline = f.readline()
    
    writefile2.write(firstline)
    writefile2.write(secondline)
    writefile3.write(firstline)
    
#     for line in f:
#         if(len(line.split())==2):
#             writefile2.write(line)
#         elif(len(line.split())==3):
#             writefile3.write(" ".join(line.split()[1:]))
#         else:
#             