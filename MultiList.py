from Tkinter import *

import General

class Table(object):

    def __init__(self, Master, nCol, ColNames, ColWidth, Spacer, Highlight, Font, Color):
        
        self.Master = Master

        # Font
        self.Font = Font
        # Font Color
        self.Color = Color

        # Scalar
        self.nCol = nCol

        self.current = None
        
        # Lists
        self.ColNames = ColNames
        self.ColWidth = ColWidth
        self.Spacer = Spacer
        self.Highlight = Highlight

        # Dictionary 
        self.Columns = dict()

    # Draws the whole Widget
    def Draw(self):

        self.build_Column()
        self.draw_Header()
        self.draw_VSB()
        self.draw_List()

    # Builds the columns (Frames)
    def build_Column(self):
        
        # Builds the frame columns
        for i in range(0, self.nCol):

            self.Columns[self.ColNames[i]] = { 'Frame': Frame(self.Master, 
                                                              width=self.ColWidth[i],
                                                              relief=RAISED,
                                                              border=1) }

            self.Columns[self.ColNames[i]]['Frame'].pack(fill=BOTH, expand=True, side=LEFT)
            self.Columns[self.ColNames[i]]['Frame'].pack_propagate(0)
            
            self.Columns[self.ColNames[i]]['Highlight'] = self.Highlight[i]
            self.Columns[self.ColNames[i]]['Spacer'] = self.Spacer[i]

    # Draws the header
    def draw_Header(self):
        
        for i in range(0, self.nCol):

            self.Columns[self.ColNames[i]]['Header'] = Label(self.Columns[self.ColNames[i]]['Frame'], 
                                                             text=self.ColNames[i], 
                                                             font=self.Font,
                                                             width=self.ColWidth[i],
                                                             relief=RAISED)

            self.Columns[self.ColNames[i]]['Header'].pack(side=TOP)
            self.Columns[self.ColNames[i]]['Header'].bind('<Button-1>', lambda event, by=self.ColNames[i]: self.Sort(event, by))

    # Draws the header
    def draw_VSB(self):

        self.vsb = Scrollbar(self.Columns[self.ColNames[self.nCol-1]]['Frame'], 
                             orient='vertical', 
                             command=self.OnVsb)

        self.vsb.pack(fill=Y, side=RIGHT)

    # Draws the Lists
    def draw_List(self):

        for i in range(0, self.nCol):

            self.Columns[self.ColNames[i]]['List'] = Listbox(self.Columns[self.ColNames[i]]['Frame'],
                                                             yscrollcommand=self.vsb.set,
                                                             selectmode=SINGLE,
                                                             selectborderwidth=0,
                                                             selectbackground=self.Color,
                                                             selectforeground='white',
                                                             highlightthickness=0,
                                                             width=self.ColWidth[i],
                                                             font=self.Font)

            self.Columns[self.ColNames[i]]['List'].pack(fill=Y, expand=True, side=LEFT)
            self.Columns[self.ColNames[i]]['List'].bind('<Button-1>', 
                                                        lambda event, List=self.Columns[self.ColNames[i]]['List']: 
                                                        self.OnButtonClick(event, List))


            # Mac/Windows
            self.Columns[self.ColNames[i]]['List'].bind('<MouseWheel>', self.OnListboxMouseWheel)

            # Linux MouseWheel Down
            self.Columns[self.ColNames[i]]['List'].bind('<Button-4>', self.OnListboxMouseWheel)
            # Linux MouseWheel Up
            self.Columns[self.ColNames[i]]['List'].bind('<Button-5>', self.OnListboxMouseWheel)

            self.Columns[self.ColNames[i]]['StringVar'] = StringVar()
            self.Columns[self.ColNames[i]]['StringVar'].set('')

    ''' ==================================================================================
    FUNCTION OnVsb: Permit to Scroll listboxes at the same time
    ==================================================================================  '''            
    def OnVsb(self, *args):

        for col in self.Columns.keys():
            self.Columns[col]['List'].yview(*args)

    ''' ==================================================================================
    FUNCTION OnButtonClick: Selects identical index from the other lists
    ==================================================================================  '''            
    def OnButtonClick(self, event, List):
        
        try:
            Index = List.nearest(event.y)
        except:
            return
        
        if Index != self.current and Index != '':
            
            if self.current != None:
                for col in self.Columns.keys():
                    if self.Columns[col]['Highlight']:
                        self.Columns[col]['List'].itemconfig(self.current, bg='white', foreground='black')
                    
            for col in self.Columns.keys():
                if self.Columns[col]['Highlight']:
                    self.Columns[col]['List'].itemconfig(Index, bg=self.Color, foreground='white')
                
                self.Columns[col]['StringVar'].set(self.Columns[col]['List'].get(Index)[1:])
        
            self.current = Index
            
    ''' ==================================================================================
    FUNCTION OnListboxMouseWheel: Scroll the Listboxes based on the mouse wheel event.
    ==================================================================================  '''
    def OnListboxMouseWheel(self, event):

        # Convert mousewheel motion to scrollbar motion.
        
        if event.num == 4:    # Linux encodes wheel as 'buttons' 4 and 5
            delta = -1
        elif event.num == 5:
            delta = 1
        else:                   # Windows & OSX
            delta = event.delta
            
        for col in self.Columns:
            self.Columns[col]['List'].yview("scroll", delta, "units")
            
        # Return 'break' to prevent the default bindings from
        # firing, which would end up scrolling the widget twice.
        return "break"

    ''' ==================================================================================
    FUNCTION Clear: Clears all the data in the table
    ==================================================================================  '''            
    def Clear(self):

        for col in self.Columns.keys():
	    try:
                self.Columns[col]['List'].selection_clear(0, END)
	    except:
		continue

        for col in self.Columns.keys():
	    try:
                self.Columns[col]['List'].delete(0, END)
	    except:
		continue
        
        self.current = None
        
    ''' ==================================================================================
    FUNCTION Add: Adds ONE item to listboxes
    ==================================================================================  '''            
    def Add(self, Item, BGColor):
        
        for i in range(0, self.nCol):
            try:
                self.Columns[self.ColNames[i]]['List'].insert(END, General.repeat(' ', self.Columns[self.ColNames[i]]['Spacer']) + str(Item[i]))
                
                if BGColor[i] != None:
                    self.Columns[self.ColNames[i]]['List'].itemconfig(self.Columns[self.ColNames[i]]['List'].size()-1, bg=BGColor[i])
            except:
                pass
            
    ''' ==================================================================================
    FUNCTION Add_List: Adds MULTIPLE items to listboxes
    ==================================================================================  '''            
    def Add_List(self, Items, BGColors):
        
        for n in range(0, len(Items)):
            for i in range(0, self.nCol):
                self.Columns[self.ColNames[i]]['List'].insert(END, General.repeat(' ', self.Columns[self.ColNames[i]]['Spacer']) + str(Items[n][i]))

                if BGColors[n][i] != None:
                    self.Columns[self.ColNames[i]]['List'].itemconfig(self.Columns[self.ColNames[i]]['List'].size()-1,
                                                                      bg=BGColors[n][i])

    ''' ==================================================================================
    FUNCTION Delete: Deletes ONE item of the listboxes
    ==================================================================================  '''            
    def Delete(self, Item, Col):
        
        Found = False
        # Does the column exist
        for i in range(0, self.nCol):
            if self.ColNames[i] == Col:
                Found = True
                break

        if Found:
            # Does item exist in column
            Index = 0
            Found = False
            for item in self.Columns[Col]['List'].get(0, END):
                if item.lstrip() == str(Item):
                    Found = True
                    break
                Index += 1

            if Found:
                # Delete all item at index in each column
                for i in range(0, self.nCol):
                    self.Columns[self.ColNames[i]]['List'].delete(Index)

    ''' ==================================================================================
    FUNCTION Set: Sets ONE item of the listboxes
    ==================================================================================  '''            
    def Set(self, Item, Col, Value, UptCol):
        
        Found = False
        # Does the column exist
        for i in range(0, self.nCol):
            if self.ColNames[i] == Col:
                Found = True
                break
            
        if Found:
            # Does item exist in column
            Index = 0
            Found = False
            for item in self.Columns[Col]['List'].get(0, END):
                #if item.replace(General.repeat(' ',self.Columns[Col]['Spacer']),'') == str(Item):
                if item.lstrip() == str(Item):
                    Found = True
                    break
                Index += 1

            if Found:
                # Does the update column exist
                Found = False
                for i in range(0, self.nCol):
                    if self.ColNames[i] == UptCol:
                        Found = True
                        break
                
                if Found:
                    # Set item value
                    self.Columns[UptCol]['List'].delete(Index)
                    
                    try:
                        self.Columns[UptCol]['List'].insert(Index, General.repeat(' ',self.Columns[UptCol]['Spacer']) + str(Value))
                    except:
                        self.Columns[UptCol]['List'].insert(END, General.repeat(' ',self.Columns[UptCol]['Spacer']) + str(Value))

    ''' ==================================================================================
    FUNCTION Sort: Sorts by a column name
    ==================================================================================  '''            
    def Sort(self, event, by):
        
        return
