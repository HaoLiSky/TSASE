import sys
import gtk
import pango
import code

yellow  ="#b58900"
orange  ="#cb4b16"
red     ="#dc322f"
magenta ="#d33682"
violet  ="#6c71c4"
blue    ="#268bd2"
cyan    ="#2aa198"
green   ="#859900"
gray    ="#657b83"
offwhite="#FDF6E3"


class ColorWriter:
    def __init__(self, textbuffer, color):
        self.textbuffer = textbuffer
        self.tag = self.textbuffer.create_tag()
        self.tag.set_property("foreground", color)
    def write(self, text):
        end = self.textbuffer.get_end_iter()
        self.textbuffer.insert_with_tags(end, text, self.tag)
        

class Console(gtk.ScrolledWindow):

    def __init__(self, callback = None, locals = None, banner = ""):
        gtk.ScrolledWindow.__init__(self)
        self.interpreter = code.InteractiveInterpreter(locals)
        self.resetbuffer()
        self.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.textview = gtk.TextView()
        self.textview.connect("key_press_event", self.event_key_pressed)
        self.textview.modify_font(pango.FontDescription("monospace 9"))
        self.textview.modify_base(0, gtk.gdk.color_parse(offwhite))
        self.textview.modify_text(0, gtk.gdk.color_parse(blue))
        self.textview.set_wrap_mode(gtk.WRAP_WORD_CHAR)
        self.textview.set_left_margin(4)
        self.textview.connect("event-after", self.event_after)
        self.textbuffer = self.textview.get_buffer()
        self.errorWriter = ColorWriter(self.textbuffer, magenta)
        self.inputWriter = ColorWriter(self.textbuffer, blue)
        self.outputWriter = ColorWriter(self.textbuffer, gray)
        self.outputWriter.write(banner)
        self.inputWriter.write(">>> ")
        self.scroll_to_bottom()
        self.promptoffset = self.textbuffer.get_end_iter().get_offset()
        self.textbuffer.connect("mark-set", self.event_mark_set)
        self.textbuffer.connect("changed", lambda w: self.textview.queue_draw)
        self.add(self.textview)
        self.history = [""]
        self.history_index = 0
        self.callback = callback

    def resetbuffer(self):
        self.buffer = []

    def push(self, line):
        self.buffer.append(line)
        source = "\n".join(self.buffer)
        sys.stdout = self.outputWriter
        sys.stderr = self.errorWriter
        more = self.interpreter.runsource(source, "<console>")
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self.scroll_to_bottom()
        if not more:
            self.resetbuffer()
            if self.callback is not None:
                self.callback()
        return more

    def scroll_to_bottom(self):
        if self.get_property("visible"):
            self.textview.scroll_to_mark(self.textbuffer.get_insert(), 0.0)

    def event_mark_set(self, tb, iter, textmark):
        if iter.get_offset() < self.promptoffset:
            self.textbuffer.move_mark(textmark, self.textbuffer.get_iter_at_offset(self.promptoffset))

    def get_item(self, item):
        if item in self.interpreter.locals:
            return self.interpreter.locals[item]
        else:
            return None
    
    def set_item(self, name, item):
        self.interpreter.locals[name] = item

    def get_current_command(self):
        start = self.textbuffer.get_iter_at_offset(self.promptoffset)
        end = self.textbuffer.get_end_iter()
        source = self.textbuffer.get_text(start, end, False)
        return source
        
    def delete_current_command(self):
        start = self.textbuffer.get_iter_at_offset(self.promptoffset)
        end = self.textbuffer.get_end_iter()
        self.textbuffer.delete(start, end)
    

    def event_key_pressed(self, widget, event):
        key = gtk.gdk.keyval_name(event.keyval)
        if key == "Up":
            self.history_index = max(0, self.history_index - 1)
            self.delete_current_command()
            self.inputWriter.write(self.history[self.history_index])
            self.scroll_to_bottom()
            return True
        if key == "Down":
            self.history_index = min(len(self.history) - 1, self.history_index + 1)
            self.delete_current_command()
            self.inputWriter.write(self.history[self.history_index])
            self.scroll_to_bottom()
            return True
        if key == "Return":
            source = self.get_current_command()
            self.history[len(self.history)-1] = source
            self.history.append("")
            self.history_index = len(self.history) - 1
            self.textbuffer.insert(self.textbuffer.get_end_iter(), "\n")
            more = self.push(source)
            self.textview.queue_draw()
            if more:
                self.inputWriter.write("... ")
            else:
                self.inputWriter.write(">>> ")
            self.promptoffset = self.textbuffer.get_end_iter().get_offset()
            self.scroll_to_bottom()
            return True
        
    
    def event_after(self, widget, event):
        if event.type == gtk.gdk.KEY_PRESS:
            key = gtk.gdk.keyval_name(event.keyval)
            if key in ["Up", "Down", "Return"]:
                return
            self.history_index = len(self.history) - 1
            self.history[self.history_index] = self.get_current_command()     
            

if __name__ == "__main__":
    w = gtk.Window(gtk.WINDOW_TOPLEVEL)
    c = Console()    
    w.add(c)
    w.show_all()
    gtk.main()

        


