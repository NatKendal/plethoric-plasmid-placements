import sys

def readingFiles(filename):
    file = open(filename, "r")
    print(file.read())
    file.close()

if __name__ == "__main__":
    print(sys.argv)
    readingFiles(sys.argv[1])
