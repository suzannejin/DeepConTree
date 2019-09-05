




def __rr_to_contactbench(prediction,native,output):
    
    
    natives = []
    with open(native) as f:
        for line in f:
            line = line.strip().split(" ")
            i = int(line[0])
            j = int(line[1])
            natives.append( (i,j) )
    
    
    predictions = []
    with open(prediction) as f:
        for line in f:
            line = line.strip().split(" ")
            i = int(line[0])
            j = int(line[1])
            p = float(line[4])
            if (i,j) in natives or (j,i) in natives:
                label = "TRUE"
            else:
                label = "FALSE"
            predictions.append( (i,j,p,label) )
            
            
    out = open(output,"wt")
    for v in sorted(predictions,key=lambda x:x[2],reverse=True):
        out.write("{} {} {} {} {}\n".format(v[0],v[1],v[1]-v[0],v[2],v[3]))
    out.close()
    
    
    
