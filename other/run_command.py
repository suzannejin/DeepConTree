



def __run_command(description_prog,command,separator):
    ''' Print and run command. '''
    
    import os
    import datetime
    
    print(separator)
    print("# {}\n# {}".format(description_prog,datetime.datetime.today().strftime('%b %d %Y %H:%M:%S')))
    print(command)
    print(" ")
    os.system(command)
