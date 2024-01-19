
def checker(num):
    
    out = ''
    
    if num %3!=0 and num%5!=0:
        return num

    if num%3==0:
        out+='fizz'
        
    if num%5==0:
        out+='buzz'
        
    return out
            
    


def main():
    
    maximum = int(input('enter max :'))
    
    for i in range(1, maximum+1):
        print(checker(i))

main()
        