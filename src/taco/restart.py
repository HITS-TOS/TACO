class Taco:
    """Test class with restart functionality"""

    def __init__(self):
        self.data = {}

    def task(self):
        self.data["b"] = 2


taco = Taco()
print(taco.data)

taco.task()
print(taco.data)
