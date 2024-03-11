





print("Stuff in global scope like this gets called every time!")

def quickFunction(s):
    print("Stuff in here only gets called if the function is called!")
    print(s)

if __name__ == "__main__":
    print("Stuff in a 'if __name__ == \"__main__\":' block only gets called if it's the main file!")

    quickFunction("And this is where you should call your functions!")
