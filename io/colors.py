
def read_colors(filename):
    lines = open(filename, 'r').readlines()
    colors = []
    index = 0
    while True:
        try:
            frameColors = []
            count = int(lines[index].strip())
            index += 1
            for i in range(count):
                color = [float(c) for c in lines[index].strip().split()]
                if len(color) < 3:
                    raise
                frameColors.append(color)
                index += 1
            colors.append(frameColors)
        except:
            if len(colors) == 0:
                raise IOError, "Could not read con file." + "\n" + str(e)
            break
    return colors
    
if __name__ == "__main__":
    print read_colors('test.colors')
    
