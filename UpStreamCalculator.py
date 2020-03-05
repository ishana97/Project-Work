from tkinter import *
from PIL import Image, ImageTk

myWindow = Tk()
myWindow.title('UpStream Cash Flow Calculator')
myWindow.geometry('1920x1080')  # See if you can do this by percentage of screen (reactive scaling)


######## SET ALL VARIABLES HERE ########

# Initialize DoubleVars
houseVar = DoubleVar(); pickupVar = DoubleVar(); trucksVar = DoubleVar(); deliveryVar = DoubleVar(); recycleVar = DoubleVar()
marketplaceVar = DoubleVar(); avgmilesdayVar = DoubleVar(); dieselcostVar = DoubleVar(); mpgVar = DoubleVar()
numemployeesVar = DoubleVar(); wageVar = DoubleVar(); timeperVar = DoubleVar(); buffertimeVar = DoubleVar();
dumpCost = DoubleVar(); housetrash = DoubleVar()

# Set Default Values
houseVar.set("80.0"); pickupVar.set("17.0"); trucksVar.set("1.0"); deliveryVar.set("2.40")
recycleVar.set("1.08"); marketplaceVar.set("1.10"); avgmilesdayVar.set("50"); dieselcostVar.set("2.90")
mpgVar.set("3.5"); numemployeesVar.set("0.333"); wageVar.set("18.00"); timeperVar.set("30")
buffertimeVar.set("2.0"); dumpCost.set("68");housetrash.set("11.4")

######## INITIALIZE ALL FRAMES HERE ########
# Left Frame Initialize
leftFrame = Frame(myWindow, width=200, height = 600)
leftFrame.grid(row=0, column=0, padx=10, pady=2)

# MIDDLE FRAME MAIN INITIALIZE
midFrameMain= Frame(myWindow, width=200, height = 1000, bg="#50BE95")
midFrameMain.grid(row=0, column=1, padx=10, pady=2)

# Mid Frame Top Initialize
midFrameTop= Frame(midFrameMain, width=200, height = 600,borderwidth=1,relief=SUNKEN)
midFrameTop.grid(row=0, column=0, padx=10, pady=70)

# Mid Frame Mid Initialize
midFrameMid= Frame(midFrameMain, width=200, height = 600,borderwidth=1,relief=SUNKEN)
midFrameMid.grid(row=1, column=0, padx=10, pady=70)

# Mid Frame Bottom Initialize
midFrameBottom= Frame(midFrameMain, width=200, height = 600,borderwidth=1,relief=SUNKEN)
midFrameBottom.grid(row=2, column=0, padx=10, pady=70)

# RIGHT FRAME MAIN INITIALIZE
rightFrameMain = Frame(myWindow, width=200, height = 1000)
rightFrameMain.grid(row=0, column=2, padx=10, pady=2)
Label(rightFrameMain, text="Positive Cash Flow",font='Helvetica 18 bold').grid(row=0, column=0,pady=2)
Label(rightFrameMain, text="Negative Cash Flow",font='Helvetica 18 bold').grid(row=2, column=0)

# Right Frame Top Initialize
rightFrameTop = Frame(rightFrameMain, width=200, height = 600)
rightFrameTop.grid(row=1, column=0, padx=10, pady=2)

# Right Frame Bottom Initialize
rightFrameBottom = Frame(rightFrameMain, width=200, height = 600)
rightFrameBottom.grid(row=3, column=0, padx=10, pady=2)

# Right Frame Bottom Bottom Initialize
rightFrameBottom2 = Frame(rightFrameMain, width=200, height = 600)
rightFrameBottom2.grid(row=4, column=0, padx=10, pady=2)

######## DEFINE ALL FUNCTIONS HERE  ########

def func(event):
    print("You hit return.")

def onclick(event=None):
    ##### Right Top Frame - Positive Cash Flow ########
    aa = pickupVar.get()
    bb = houseVar.get()
    pickupcf = round(float(aa) * float(bb), 2)
    Label(rightFrameTop, width=6,bg ="white", text=pickupcf).grid(row=1, column=2)

    cc = deliveryVar.get()
    delcf = round(float(cc) * float(bb), 2)
    Label(rightFrameTop, width=6, bg="white", text=delcf).grid(row=2, column=2)

    dd = recycleVar.get()
    reccf = round(float(dd) * float(bb), 2)
    Label(rightFrameTop, width=6, bg="white", text=reccf).grid(row=3, column=2)

    ee = marketplaceVar.get()
    mrkcf = round(float(ee)*float(bb),2)
    Label(rightFrameTop, width=6, bg="white", text=mrkcf).grid(row=4, column=2)

    posflow = round(pickupcf+delcf+reccf+mrkcf,2)
    Label(rightFrameTop, width=6, bg="#21946a", text=posflow).grid(row=5, column=2)

    ##### Mid Top Frame - Fuel Calculation ########
    ff = avgmilesdayVar.get(); gg = dieselcostVar.get(); hh = mpgVar.get();
    fuelcostday = round((float(ff)/float(hh))*float(gg),2)
    Label(midFrameTop, width=7, relief=GROOVE, text=fuelcostday).grid(row=4, column=1)
    #############

    if ((round((bb/300)*2)) % 2) == 0:
        routes = round(((bb/300)*2))
    else:
        routes = round(((2*(bb/300))+1))


    Label(leftFrame, width=7, relief=GROOVE, text=routes).grid(row=6, column=0)

    ######## Mid Frame Mid - Labor Calculation ########

    unitemployee = round((bb / routes) / 150, 2)
    Label(midFrameMid, width=7, bg="gray", text=unitemployee).grid(row=1, column=1)

    ii = timeperVar.get()
    timedaytotal = round(float(bb) * (float(ii) / 60 / 60), 2)
    Label(midFrameMid, width=7, bg="gray", text=timedaytotal).grid(row=5, column=1)

    jj = buffertimeVar.get()
    kk = wageVar.get()
    labordaytotal = round((timedaytotal + float(jj)) * float(kk) * unitemployee, 2)
    Label(midFrameMid, width=7, bg="white", text=labordaytotal).grid(row=6, column=1)

    ######## Right Frame Bottom - Negative Cash Flow########
    fueltotal = round(fuelcostday*routes,2)
    Label(rightFrameBottom, width=6,bg="white", text=fueltotal).grid(row=1, column=1)

    maintainence = float(500)
    Label(rightFrameBottom, width=6, bg="white", text=maintainence).grid(row=2, column=1)

    insurance = float(1250)
    Label(rightFrameBottom, width=6, bg="white", text=insurance).grid(row=3, column=1)

    labortotal = round(labordaytotal * routes,2)
    Label(rightFrameBottom, width=6, bg="white", text=labortotal).grid(row=4, column=1)

    ll = dumpCost.get(); mm = housetrash.get()
    dumptotalcost = round(((float(mm)*30)/2000)*ll*bb,2)
    Label(rightFrameBottom, width=6, bg="white", text=dumptotalcost).grid(row=5, column=1)

    negflow = round(fueltotal + maintainence + insurance + labortotal + dumptotalcost,2)
    Label(rightFrameBottom, width=6, bg="#d94e61", text=negflow).grid(row=6, column=1)

    ######## Main Cash Flow Values ########

    cashflow = round (posflow - negflow,2)
    Label(rightFrameBottom2, width=10, bg="white", text=cashflow,font='Helvetica 18 bold').grid(row=0, column=0,padx=10)

myWindow.bind('<Return>', onclick) # Enables using 'enter' button to update fields

######## Left Frame Contents########

# Logo
imagename="Capture_burned.png"
im1 = Image.open(imagename)
size = (im1.width//10,im1.height//10)
im1 = ImageTk.PhotoImage(Image.open(imagename).resize(size))
Label(leftFrame,image=im1).grid(row=0,column=0)

# Labels & Spinbox
Label(leftFrame, text="# Houses").grid(row=1, column=0)
Label(leftFrame, text="# Trucks").grid(row=3, column=0)
Label(leftFrame, text="Routes Per Week").grid(row=5, column=0)

Label(leftFrame, width=7, relief=GROOVE, text="2").grid(row=6, column=0)

Spinbox(leftFrame, from_=0, to=2000, increment=1, width=5, textvariable=houseVar, command=onclick).grid(row=2, column=0)
Spinbox(leftFrame, from_=0, to=100, increment=1, width=5, textvariable=trucksVar).grid(row=4, column=0)


######## Middle Frame Top Contents ########

# Positive Cash Flow Section
Label(rightFrameTop, text="Waste Collection").grid(row=1, column=0)
Label(rightFrameTop, text="Order Delivery").grid(row=2, column=0)
Label(rightFrameTop, text="Selling Recyclables").grid(row=3, column=0)
Label(rightFrameTop, text="Second Hand Marketplace").grid(row=4, column=0)
Label(rightFrameTop, text="Per Household").grid(row=0, column=1)

Label(rightFrameTop, width=6, bg="gray", text="17").grid(row=1, column=2)
Label(rightFrameTop, width=6, bg="gray", text="2.40").grid(row=2, column=2)
Label(rightFrameTop, width=6, bg="gray", text="1.08").grid(row=3, column=2)
Label(rightFrameTop, width=6, bg="gray", text="1.10").grid(row=4, column=2)
Label(rightFrameTop, width=6, bg="green", text="1726.4").grid(row=5, column=2)

Spinbox(rightFrameTop, from_=0, to=100, increment=0.1, width=5, textvariable=pickupVar, command=onclick).grid(row=1, column=1)
Spinbox(rightFrameTop, from_=0, to=100, increment=0.01, width=5, textvariable=deliveryVar, command=onclick).grid(row=2, column=1)
Spinbox(rightFrameTop, from_=0, to=100, increment=0.01, width=5, textvariable=recycleVar, command=onclick).grid(row=3, column=1)
Spinbox(rightFrameTop, from_=0, to=100, increment=0.01, width=5, textvariable=marketplaceVar, command=onclick).grid(row=4, column=1)


######## Right Frame Bottom Contents ########

# Negative Cash Flow Labels
Label(rightFrameBottom, text="Fuel Cost").grid(row=1, column=0, padx =100)
Label(rightFrameBottom, text="Maintainence").grid(row=2, column=0)
Label(rightFrameBottom, text="Insurance").grid(row=3, column=0)
Label(rightFrameBottom, text="Labor").grid(row=4, column=0)
Label(rightFrameBottom, text="Dumping Fees").grid(row=5, column=0)
Label(rightFrameBottom, text="Negative Cash Flow Total").grid(row=6, column=0)

Label(rightFrameBottom, width=6,bg="gray", text="82.86").grid(row=1, column=1,padx=50)
Label(rightFrameBottom, width=6, bg="gray", text="500").grid(row=2, column=1)
Label(rightFrameBottom, width=6, bg="gray", text="1250").grid(row=3, column=1)
Label(rightFrameBottom, width=6, bg="gray", text="25.96").grid(row=4, column=1)
Label(rightFrameBottom, width=6, bg="gray", text="930.24").grid(row=5, column=1)
Label(rightFrameBottom, width=6, bg="red", text="2789.06").grid(row=6, column=1)

######## Middle Frame Top Contents ########

# Fuel Calculation Box

Label(midFrameTop, text="Fuel Information",font='Helvetica 18 bold').grid(row=0, column=0)
Label(midFrameTop, text="Average Miles Driven Per Day").grid(row=1, column=0)
Label(midFrameTop, text="Diesel Cost Per Gallon").grid(row=2, column=0)
Label(midFrameTop, text="Miles Per Gallon").grid(row=3, column=0)
Label(midFrameTop, text="Fuel Cost Per Day").grid(row=4, column=0)

Label(midFrameTop, width=7, relief=GROOVE, text="41.43").grid(row=4, column=1)

Spinbox(midFrameTop, from_=0, to=100, increment=1, width=5, textvariable=avgmilesdayVar, command=onclick).grid(row=1, column=1)
Spinbox(midFrameTop, from_=0, to=100, increment=0.01, width=5, textvariable=dieselcostVar, command=onclick).grid(row=2, column=1)
Spinbox(midFrameTop, from_=0, to=100, increment=0.01, width=5, textvariable=mpgVar, command=onclick).grid(row=3, column=1)


######## Middle Frame Mid Contents ########

# Labor Calculation Box

Label(midFrameMid, text="Labor Information",font='Helvetica 18 bold').grid(row=0, column=0)
Label(midFrameMid, text="Unit Employees").grid(row=1, column=0)
Label(midFrameMid, text="Hourly Rate").grid(row=2, column=0)
Label(midFrameMid, text="Time per House [s]").grid(row=3, column=0)
Label(midFrameMid, text="Buffer Time [Hrs]").grid(row=4, column=0)
Label(midFrameMid, text="Total Daily Time [Hrs]").grid(row=5, column=0)
Label(midFrameMid, text="Labor Cost per Day").grid(row=6, column=0)

Label(midFrameMid, width=7, bg="gray", text="0.27").grid(row=1, column=1)
Label(midFrameMid, width=7, bg="gray", text="0.67").grid(row=5, column=1)
Label(midFrameMid, width=7, bg="gray", text="12.98").grid(row=6, column=1)

Spinbox(midFrameMid, from_=0, to=100, increment=0.1, width=5, textvariable=wageVar, command=onclick).grid(row=2, column=1)
Spinbox(midFrameMid, from_=0, to=100, increment=1.0, width=5, textvariable=timeperVar, command=onclick).grid(row=3, column=1)
Spinbox(midFrameMid, from_=0, to=10, increment=0.5, width=5, textvariable=buffertimeVar, command=onclick).grid(row=4, column=1)

######## Middle Frame Bottom Contents ########

# Dumping Fees Box

Label(midFrameBottom, text="Dumping Information", font='Helvetica 18 bold').grid(row=0, column=0)
Label(midFrameBottom, text="County Dumping Costs[$/tons]").grid(row=1, column=0)
Label(midFrameBottom, text="Avg Household Trash Per Day [lbs]").grid(row=2, column=0)

Spinbox(midFrameBottom, from_=0, to=100, increment=0.1, width=5, textvariable=dumpCost, command=onclick).grid(row=1, column=1)
Spinbox(midFrameBottom, from_=0, to=100, increment=0.1, width=5, textvariable=housetrash, command=onclick).grid(row=2, column=1)

myWindow.mainloop()









