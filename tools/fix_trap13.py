import csv

# trap 13 has an extra column
# run this once to remove said column
# will need to give it a new file to output to, then can delete the original file and name it trap13 again :)

def deleteColumn17(filename,tempName):

    with open(filename,"r") as file:

        reader = csv.reader(file)
    
        with open(tempName,"w",newline='') as fixed:

            writer = csv.writer(fixed)
            
            for row in reader:
                
                writer.writerow(tuple(row[i] for i in range(len(row)) if i != 17))