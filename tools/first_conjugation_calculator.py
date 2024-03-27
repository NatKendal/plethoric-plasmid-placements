import csv

# ------------------------------------------------------------------
### Checking the time step of the first conjugation event in each file
# ------------------------------------------------------------------

# input: file name
# output: time step of first detected conjugation event
def getFirstConj(filename):
    
    #read in the csv file and store the information in the dictionaries 
    with open(filename, 'r') as file:
        
        reader = csv.reader(file, delimiter=',')

        # since there is a header, we skip line one
        next(reader)

        for row in reader:

            # check the transconjugant flag
            # if it is one, return the current time step
            flag = int(row[21])
            if flag == 1:
                return int(row[0])
            

if __name__ == "__main__":
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap5.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap6.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap8.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap9.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap10.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap11.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap12.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap13.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap14.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap16.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap17.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap18.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap19.csv"))
    print(getFirstConj(r"C:\Users\nat\Documents\GitHub\plethoric-plasmid-placements\raw_data\old_dataset_1\trap20.csv"))
