import datetime
import csv

def write(filename,settings,*args,**kwargs):

    assert(type(filename) == str)
    assert(type(settings) == dict)

    content = open(filename+'.txt','w')

    t  = datetime.datetime.now()
    content.write('Written on:\n' + t.strftime('%x  %X\n\n'))
    content.write('______________________________\n\n')
    content.write('(Settings:)\n\n'+str(settings)+'\n\n')
    content.write('______________________________\n\n')
    content.write('(Values:)\n\n')

    for key, value in kwargs.items():
        content.write('%%'+str(key)+'%%\n')
        content.write('$$'+str(value)+'$$\n\n')
    
    content.write('\n\n\n\n')
    content.close()

    print('File ',filename+'.txt is written')

    return      

def read(filename):

    content = open(filename+'.txt','r')
    text = content.read()

    cursor   = text.find('(Values:)')

    assert text.find('(Values:)',cursor+1) <= 0, '(Values:) cannot be used as an input variable name or value'

    var = {}

    while True:

        cursor = text.find('%%',cursor+1)
        begin  = cursor+2
        cursor = text.find('%%',cursor+1)
        end    = cursor
        key    = text[begin:end]
        cursor = text.find('$$',cursor+1)
        begin  = cursor+2
        cursor = text.find('$$',cursor+1)
        end    = cursor
        val    = text[begin:end]
       
        val    = [eval(x) for x in val.strip("[").strip("]").split(",")]  

        var[key] = val

        if text.find('%%',cursor) <= 0:
            break

    return var

def appendcsv(filename,data,*args,**kwargs):

    assert(type(filename) == str)

    with open(filename, 'a') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow(data)

    csvFile.close()

    print('CSV file ',filename+'.txt is edited')

    return      
