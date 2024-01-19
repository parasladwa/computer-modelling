
import random

    
def flip_coin():
    # Generate a random number between 0 and 1
    result = random.randint(0, 1)

    # Return "Heads" for 0 and "Tails" for 1
    if result == 0:
        return "Heads"
    else:
        return "Tails"


def game():
    
    counter = 0
    winnings = 0
    
    
    coin = flip_coin()
    while coin!='Tails':
        counter+=1
        winnings+= 2**counter
        coin = flip_coin()
    return winnings

total =0
for i in range(0, 1000000):
    total += game()
print( total/1000000)
    



























