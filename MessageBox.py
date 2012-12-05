from Tkinter import *

class MessageItem(Frame):

    def __init__(self, master, message, **kwds):
        Frame.__init__(self, master, **kwds)
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.text = Label(self, text=message, anchor='w', bg='gold')
        self.text.grid(row=0, column=0, sticky='nsew')

class scrollableContainer(Frame):

    def __init__(self, master, **kwargs):
        Frame.__init__(self, master, **kwargs) #holds canvas & scrollbars
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        self.canv = Canvas(self, bd=0, highlightthickness=0)
        self.hScroll = Scrollbar(self, orient='horizontal',
                                 command=self.canv.xview)
        self.hScroll.grid(row=1, column=0, sticky='we')
        self.vScroll = Scrollbar(self, orient='vertical',
                                 command=self.canv.yview)
        self.vScroll.grid(row=0, column=1, sticky='ns')
        self.canv.grid(row=0, column=0, sticky='nsew')        
        self.canv.configure(xscrollcommand=self.hScroll.set,
                            yscrollcommand=self.vScroll.set)

        self.frm = Frame(self.canv, bd=2, bg='green') #holds messages
        self.frm.grid_columnconfigure(0, weight=1)

        self.canv.create_window(0, 0, window=self.frm, anchor='nw', tags='inner')

        self.messages = []
        for i in range(20):
            m = MessageItem(self.frm, 'Something Profound', bd=2, bg='black')
            m.grid(row=i, column=0, sticky='nsew', padx=2, pady=2)
            self.messages.append(m)

        self.update_layout()        
        self.canv.bind('<Configure>', self.on_configure)

    def update_layout(self):
        self.frm.update_idletasks()
        self.canv.configure(scrollregion=self.canv.bbox('all'))
        self.canv.yview('moveto','1.0')
        self.size = self.frm.grid_size()

    def on_configure(self, event):
        w,h = event.width, event.height
        natural = self.frm.winfo_reqwidth()
        self.canv.itemconfigure('inner', width= w if w>natural else natural)
        self.canv.configure(scrollregion=self.canv.bbox('all'))

    def add_message(self, message, col_index):
        m = MessageItem(self.frm, message, bd=2, bg='red')
        m.grid(row=self.size[1], column=0, padx=2, pady=2, sticky='we')
        self.messages.append(m)
        self.update_layout()

