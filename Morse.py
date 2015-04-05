"""
. winsound.Beep(500,150)
- winsound.Beep(500,450)
delay between Beeps : 0.1sec
delay between letters : 0.3sec
delay between words : 0.7sec
"""
import winsound
import time

def Sound(signal):
    if(type(signal) == str):
        for s in list(signal):
            if (s == '.'):
                winsound.Beep(500,150)
                time.sleep(0.1)
            elif (s == '-'):
                winsound.Beep(500,450)
                time.sleep(0.1)
        time.sleep(0.2)
    else:
        raise(TypeError)
        
"""
Switch is from:
http://stackoverflow.com/a/6606504/4706932
"""
class switch(object):
    value = None
    def __new__(class_, value):
        class_.value = value
        return True

def case(*args):
    return any((arg == switch.value for arg in args))

x = raw_input("Input the sentence:")

for char in list(x):
    
    if('qwertyuiopasdfghjklzxcvbnm 0123456789'.find(char) == -1):
        pass
    
    else:
        while switch(char):
            if case('a'):
                Sound('.-')
                break
            if case('b'):
                Sound('-...')
                break
            if case('c'):
                Sound('-.-.')
                break
            if case('d'):
                Sound('-..')
                break
            if case('e'):
                Sound('.')
                break
            if case('f'):
                Sound('..-.')
                break
            if case('g'):
                Sound('--.')
                break
            if case('h'):
                Sound('....')
                break
            if case('i'):
                Sound('..')
                break
            if case('j'):
                Sound('.---')
                break
            if case('k'):
                Sound('-.-')
                break
            if case('l'):
                Sound('..')
                break
            if case('m'):
                Sound('--')
                break
            if case('n'):
                Sound('-.')
                break
            if case('o'):
                Sound('---')
                break
            if case('p'):
                Sound('.--.')
                break
            if case('q'):
                Sound('--.-')
                break
            if case('r'):
                Sound('.-.')
                break
            if case('s'):
                Sound('...')
                break
            if case('t'):
                Sound('-')
                break
            if case('u'):
                Sound('..-')
                break
            if case('v'):
                Sound('...-')
                break
            if case('w'):
                Sound('.--')
                break
            if case('x'):
                Sound('-..-')
                break
            if case('y'):
                Sound('-.--')
                break
            if case('z'):
                Sound('--..')
                break
            if case('1'):
                Sound('.----')
                break
            if case('2'):
                Sound('..---')
                break
            if case('3'):
                Sound('...--')
                break
            if case('4'):
                Sound('....-')
                break
            if case('5'):
                Sound('.....')
                break
            if case('6'):
                Sound('-....')
                break
            if case('7'):
                Sound('--...')
                break
            if case('8'):
                Sound('---..')
                break
            if case('9'):
                Sound('----.')
                break
            if case('0'):
                Sound('-----')
                break
            if case(' '):
                time.sleep(0.7)
