class Car:
    def __init__(self,color,speed,type):
        self.color = color
        self.speed = speed
        self.type = type

    def getColor(self):
        return self.color

    def calculateTime(self):
        distance = int(input("Enter the distance: "))
        time =  distance/self.speed 
        return time

audi = Car("green",20,"Sport")
print(audi.getColor())
print(audi.calculateTime(),"hours")
