#!/usr/bin/env python

import os
import sys
import math
import time

import gtk
import gtk.glade
import gtk.gdk as gdk
import gobject
import pango
import cairo

import numpy as np

import ase
from ase.io import vasp
import tsase
from tsase.data import *
from console import Console


class queueitem:
    def __init__(self, kind):
        self.kind = kind


banner = \
"""
Python console for tsase-xyz
----------------------------
The variable p in this console is the ase.atoms.Atoms object or list of them that is being displayed. The modules ase and tsase have already been imported.
"""

class xyz(gtk.Window):

    def __init__(self):
        gtk.Window.__init__(self, gtk.WINDOW_TOPLEVEL)
        self.acquire_widgets()
        self.initialize_console()
        self.connect_events()
        self.add(self.gladewindow)
        self.set_resizable(True)
        self.moviescale.set_increments(1, 1)
        self.statusbar.modify_font(pango.FontDescription("monospace 10"))
        self.area.set_size_request(128, 128)
        self.set_default_size(640, 512)
        self.show()
        self.initialize_members()
        self.event_configure()
        self.gfx_setup_colors()
        self.gfx_reset_transform()
        gobject.timeout_add(8, self.event_timeout)

    def initialize_console(self):
        self.console = Console(callback=self.event_console_command, banner = banner)
        self.console.push("from ase import *\n")
        self.console.push("import tsase\n")
        self.consolebox.add(self.console)

    def acquire_widgets(self):
        filename = os.path.join(os.path.dirname(__file__), "xyz.glade")
        gladetree = gtk.glade.XML(filename)
        self.gladewindow           = gladetree.get_widget("window")
        self.moviescale            = gladetree.get_widget("moviescale")
        self.moviebox              = gladetree.get_widget("moviebox")
        self.playbutton            = gladetree.get_widget("playbutton")
        self.playimage             = gladetree.get_widget("playimage")
        self.boxbutton             = gladetree.get_widget("boxbutton")
        self.axisbutton            = gladetree.get_widget("axisbutton")
        self.resetbutton           = gladetree.get_widget("resetbutton")
        self.zoombutton            = gladetree.get_widget("zoombutton")
        self.radiusbutton          = gladetree.get_widget("radiusbutton")
        self.frozenbutton          = gladetree.get_widget("frozenbutton")
        self.movebutton            = gladetree.get_widget("movebutton")
        self.fps                   = gladetree.get_widget("fps")
        self.repeatx               = gladetree.get_widget("repeatx")
        self.repeaty               = gladetree.get_widget("repeaty")
        self.repeatz               = gladetree.get_widget("repeatz")
        self.holder                = gladetree.get_widget("holder")
        self.consolebox            = gladetree.get_widget("consolebox")
        self.statusbar             = gladetree.get_widget("statusbar")
        self.area                  = gladetree.get_widget("atomview")
        self.menuFileOpen          = gladetree.get_widget("menuFileOpen")
        self.menuFileOpenView      = gladetree.get_widget("menuFileOpenView")
        self.menuFileSaveAs        = gladetree.get_widget("menuFileSaveAs")
        self.menuFileSaveView      = gladetree.get_widget("menuFileSaveView")
        self.menuFileExport        = gladetree.get_widget("menuFileExport")
        self.menuFileQuit          = gladetree.get_widget("menuFileQuit")
        self.menuToolsConsole      = gladetree.get_widget("menuToolsConsole")
        self.menuHelpDocumentation = gladetree.get_widget("menuHelpDocumentation")

    def connect_events(self):
        self.connect("destroy", self.event_close)
        self.connect("key_release_event", self.event_key_released)
        self.connect("key_press_event", self.event_key_pressed)
        self.area.connect("expose_event", self.event_exposed)
        self.area.connect("configure_event", self.event_configure)
        self.area.connect("button_press_event", self.event_button_press)
        self.area.connect("button_release_event", self.event_button_release)
        self.area.connect("motion_notify_event", self.event_mouse_move)
        self.area.connect("scroll_event", self.event_scroll)
        self.repeatx.connect("value-changed", lambda w: self.gfx_center_atoms())
        self.repeaty.connect("value-changed", lambda w: self.gfx_center_atoms())
        self.repeatz.connect("value-changed", lambda w: self.gfx_center_atoms())
        self.playbutton.connect("clicked", self.event_toggle_play)
        self.resetbutton.connect("clicked", lambda w: self.gfx_reset_transform())
        self.boxbutton.connect("clicked", lambda w: self.queue_draw())
        self.axisbutton.connect("clicked", lambda w: self.queue_draw())
        self.frozenbutton.connect("clicked", lambda w: self.queue_draw())
        self.zoombutton.connect("value-changed", lambda w: self.queue_draw())
        self.radiusbutton.connect("value-changed", lambda w: self.queue_draw())
        self.moviescale.connect("value-changed", lambda w: self.queue_draw())
        self.menuFileOpen.connect("activate", self.event_menuFileOpen)
        self.menuFileOpenView.connect("activate", self.event_menuFileOpenView)
        self.menuFileSaveAs.connect("activate", self.event_menuFileSaveAs)
        self.menuFileSaveView.connect("activate", self.event_menuFileSaveView)
        self.menuFileExport.connect("activate", self.event_menuFileExport)
        self.menuFileQuit.connect("activate", self.event_close)
        self.menuToolsConsole.connect("activate", self.event_menuToolsConsole)
        self.menuHelpDocumentation.connect("activate", self.display_docs)

    def display_docs(self, *args):
        gladetree = gtk.glade.XML(os.path.join(os.path.dirname(__file__), "xyz.glade"))
        helpwindow = gladetree.get_widget("helpwindow")
        helpview = gladetree.get_widget("helpview")
        f = open(os.path.join(os.path.dirname(__file__), "xyz.help"), 'r')
        lines = ''.join(f.readlines())
        f.close()
        buffer = gtk.TextBuffer()
        buffer.set_text(lines)
        helpview.set_buffer(buffer)
        font = pango.FontDescription("monospace 10")
        helpview.modify_font(font)
        helpwindow.show_all()

    def initialize_members(self):
        self.last_draw = 0.0
        self.playing = False
        self.queue = []
        self.trajectory = None
        self.button1 = False
        self.button2 = False
        self.button3 = False
        self.mouselast = [None, None]
        self.background = [1.0, 1.0, 1.0]
        self.keys = {}
        self.black_gc = self.area.get_style().black_gc
        self.white_gc = self.area.get_style().white_gc
        self.background_gc = self.gfx_get_color_gc(self.background[0], self.background[1], self.background[2])
        self.pixmap = None
        self.repeat = (1, 1, 1)
        self.lastTime = time.time()
        self.screenatoms = []
        self.grabbedatom = None
        self.colors = []
        self.render_cairo = False
        self.pwd = os.getcwd()

    def gui_key_on(self, key):
        return key in self.keys
        
#
# EVENT -----------------------------------------------------------------------------------------
#


    def event_configure(self, *args):
        """ Called when the window is resized and created. """
        x, y, self.width, self.height = self.area.get_allocation()
        self.pixmap = gdk.Pixmap(self.area.window, self.width, self.height)     
        return True


    def event_key_pressed(self, widget, event):
        key = gdk.keyval_name(event.keyval)
        self.keys[key] = True
        if key == "space" and self.trajectory is not None and len(self.trajectory) > 1:
            self.event_toggle_play()
        if key == 'c' and (event.state & gdk.CONTROL_MASK):
            sys.exit()

    def event_key_released(self, widget, event):
        key = gdk.keyval_name(event.keyval)
        if self.gui_key_on(key):
            self.keys.pop(key)


    def event_mouse_move(self, widget, event):
        mx, my, mask = self.area.window.get_pointer()
        if self.mouselast[0] == None:
            self.mouselast = (mx, my)
        else:
            dx = mx - self.mouselast[0]
            dy = my - self.mouselast[1]
            self.mouselast = (mx, my)
        if self.button1:
            if self.gui_key_on("m") or (self.movebutton.get_active() and self.grabbedatom is not None):
                if self.grabbedatom is not None:
                    r = self.get_frame_atoms().get_positions() 
                    r[self.grabbedatom] += np.dot(np.linalg.inv(self.rotation), 
                                                  np.array([dx, -dy, 0]) / \
                                                  self.zoombutton.get_value() / 2.0)
                    self.get_frame_atoms().positions = r
            else:
                self.gfx_rot_x(dy * 0.009)
                self.gfx_rot_y(dx * 0.009)
            self.gfx_render()
        elif self.button2:
            self.gfx_rot_z(-dx * 0.0078125)
            self.gfx_render()
        elif self.button3:
            self.translate += np.array([dx, -dy, 0]) / self.zoombutton.get_value() / 2
            self.gfx_render()
        atomid = self.get_mouse_atom()
        if atomid is not None:     
            atom = self.get_frame_atoms()[atomid]
            r = atom.get_position()
            self.statusbar.set_text("Atom %d, %s (%.3fx %.3fy %.3fz)" % 
                                     (atomid, atom.symbol, r[0], r[1], r[2]))
        else:
            self.statusbar.set_text("")
        return True
    

    def get_mouse_atom(self):
        mx, my = self.mouselast
        atomid = None
        for a in self.screenatoms:
            d2 = (a[0] - mx)**2 + (a[1] - my)**2
            if d2 < a[2]**2:
                atomid = a[3]
        return atomid


    def event_button_press(self, widget, event):
        if event.button == 1:
            self.button1 = True
            self.grabbedatom = self.get_mouse_atom()
        if event.button == 2:
            self.button2 = True
        if event.button == 3:
            self.button3 = True
        return True


    def event_button_release(self, widget, event):
        if event.button == 1:
            self.button1 = False
            self.grabbedatom = None
        if event.button == 2:
            self.button2 = False
        if event.button == 3:
            self.button3 = False
        return True
            

    def event_scroll(self, widget, event):
        if self.gui_key_on("r"): 
            if event.direction == gdk.SCROLL_UP:
                self.radiusbutton.set_value(self.radiusbutton.get_value() * 1.1)
                self.queue_draw()
            elif event.direction == gdk.SCROLL_DOWN:
                self.radiusbutton.set_value(self.radiusbutton.get_value() * 0.9)
                self.queue_draw()
        else:
            if event.direction == gdk.SCROLL_UP:
                self.zoombutton.set_value(self.zoombutton.get_value() * 1.1)
            elif event.direction == gdk.SCROLL_DOWN:
                self.zoombutton.set_value(self.zoombutton.get_value() * 0.9)
            self.queue_draw()
        return True
                                                    

    def get_frame_atoms(self):
        drawpoint = self.trajectory[0]
        if len(self.trajectory) > 1:
            drawpoint = self.trajectory[int(self.moviescale.get_value())]
        return drawpoint
                                                                                    

    def event_exposed(self, *args):
        self.gfx_render()


    def event_toggle_play(self, *args):
        self.playing = not self.playing
        if self.playing:
            self.playimage.set_from_stock(gtk.STOCK_MEDIA_PAUSE, gtk.ICON_SIZE_BUTTON)
        else:
            self.playimage.set_from_stock(gtk.STOCK_MEDIA_PLAY, gtk.ICON_SIZE_BUTTON)
                        
            
    def event_timeout(self):
        if self.trajectory is None:
            return True
        if self.playing and len(self.trajectory) > 1:
            if time.time() - self.last_draw > 1.0/self.fps.get_value():
                nextframe = self.moviescale.get_value() + 1
                if nextframe > len(self.trajectory) - 1:
                    nextframe = 0
                self.moviescale.set_value(nextframe)
                self.queue_draw()        
        return True
    

    def event_close(self, *args):
        gtk.main_quit()
        

    def getFilenameOpen(self):
        buttons = (gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_OPEN,gtk.RESPONSE_OK)
        chooser = gtk.FileChooserDialog(title=None,action=gtk.FILE_CHOOSER_ACTION_OPEN,buttons=buttons)
        chooser.set_current_folder(self.pwd)
        response = chooser.run()
        filename = chooser.get_filename()
        chooser.destroy()
        return response, filename    

    def getFilenameSave(self):
        buttons = (gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_SAVE,gtk.RESPONSE_OK)
        chooser = gtk.FileChooserDialog(title=None,action=gtk.FILE_CHOOSER_ACTION_SAVE,buttons=buttons)
        chooser.set_current_folder(self.pwd)
        response = chooser.run()
        filename = chooser.get_filename()
        chooser.destroy()
        return response, filename    
        
    def event_menuFileOpen(self, *args):
        response, filename = self.getFilenameOpen()
        if response == gtk.RESPONSE_OK:
            self.pwd = os.path.dirname(filename)
            self.data_read(filename)
        return True

    def event_menuFileOpenView(self, *args):
        response, filename = self.getFilenameOpen()
        if response == gtk.RESPONSE_OK:
            self.pwd = os.path.dirname(filename)
            view = eval(open(filename, 'r').readline())
            self.radiusbutton.set_value(view['radius'])
            self.zoombutton.set_value(view['zoom'])
            self.translate = np.array(view['translate'])
            self.rotation = np.array(view['rotation'])
            self.repeatx.set_value(view['repeat'][0])
            self.repeaty.set_value(view['repeat'][1])
            self.repeatz.set_value(view['repeat'][2])
            self.gfx_render()
        return True

    def event_menuFileSaveAs(self, *args):
        if self.trajectory is None:
            return
        response, filename = self.getFilenameSave()
        if response == gtk.RESPONSE_OK:
            self.pwd = os.path.dirname(filename)
            self.data_write(filename)
        return True

    def event_menuFileSaveView(self, *args):
        response, filename = self.getFilenameSave()
        if response == gtk.RESPONSE_OK:
            self.pwd = os.path.dirname(filename)
            repeat = (int(self.repeatx.get_value()), 
                      int(self.repeaty.get_value()), 
                      int(self.repeatz.get_value()))
            view = repr({'radius':    self.radiusbutton.get_value(), 
                         'zoom':      self.zoombutton.get_value(), 
                         'rotation':  self.rotation.tolist(), 
                         'translate': self.translate.tolist(),
                         'repeat':    repeat})
            f = open(filename, 'w')
            f.write(view)
            f.close()
        return True

    def event_menuFileExport(self, *args):
        response, filename = self.getFilenameSave()
        if response == gtk.RESPONSE_OK:
            self.pwd = os.path.dirname(filename)
            if filename.endswith(".ps"):
                csurf = cairo.PSSurface(filename, self.width, self.height)
            if filename.endswith(".eps"):
                csurf = cairo.PSSurface(filename, self.width, self.height)
                try:
                    csurf.set_eps(True)
                except:
                    print "Could not enable eps."
            elif filename.endswith(".svg"):
                csurf = cairo.SVGSurface(filename, self.width, self.height)
            elif filename.endswith(".pdf"):
                csurf = cairo.PDFSurface(filename, self.width, self.height)
            self.cairo_context = cairo.Context(csurf)
            self.cairo_context.rectangle(0, 0, self.width, self.height)
            self.cairo_context.clip()
            self.render_cairo = True
            self.gfx_render()
            self.render_cairo = False
            csurf.finish()
        return True

    def event_menuToolsConsole(self, *args):
        if self.menuToolsConsole.get_active():
            self.consolebox.show_all()
        else:
            self.consolebox.hide()

    def event_console_command(self):
        pass
        p = self.console.get_item("p")
        if p is None:
            return
        if type(p) == ase.atoms.Atoms:
            self.data_set(p)
            return
        if type(p) == list:
            for i in p:
                if type(i) != ase.atoms.Atoms:
                    self.console.write("The p variable must be an ase.atoms.Atoms type or a list of them.\n")
                    return
            self.data_set(p)
            return
        self.console.write("The p variable must be an ase.atoms.Atoms type or a list of them.\n")
            
        

#
# GRAPHICS --------------------------------------------------------------------------------------
#

    def gfx_render(self):
        if self.trajectory is None:
            return
        self.gfx_clear()
        self.queue = []
        self.gfx_queue_atoms()
        if self.boxbutton.get_active():
            self.gfx_queue_box()
        self.gfx_transform_queue()
        self.gfx_sort_queue()
        self.gfx_draw_queue()
        if self.axisbutton.get_active():
            self.gfx_draw_axes()
        self.area.window.draw_drawable(self.white_gc, self.pixmap, 0, 0, 0, 0, self.width, self.height)
        self.last_draw = time.time()

    def gfx_get_color_gc(self, r, g, b):
        rgb = (int(r * 65535), int(g * 65535), int(b * 65535))
        gc = self.area.window.new_gc()
        gc.set_rgb_fg_color(gdk.Color(rgb[0], rgb[1], rgb[2]))
        return gc
        
    def gfx_setup_colors(self):
        for i in range(num_elements):
            c = elements[i]['color']
            self.colors.append(self.gfx_get_color_gc(c[0], c[1], c[2]))

    def gfx_reset_transform(self):
        self.radiusbutton.set_value(1.5)
        self.zoombutton.set_value(8.0)
        self.rotation = np.identity(3)
        self.translate = np.array([0.0, 0.0, 16.0])
        self.queue_draw()
    
    def gfx_queue_line(self, r1, r2, color, width = 1):
        line = queueitem("line")
        line.r1 = np.copy(r1)
        line.r2 = np.copy(r2)
        line.color = color
        line.depth = (r1 + r2) / 2.0
        line.width = width
        self.queue.append(line)
        
    def gfx_queue_atoms(self):
        ra = self.get_frame_atoms().repeat((int(self.repeatx.get_value()), 
                                            int(self.repeaty.get_value()), 
                                            int(self.repeatz.get_value())))
        r = ra.get_positions()
        symbols = ra.get_chemical_symbols()
        for i in range(len(r)):
            atom = queueitem("atom")
            atom.r = np.copy(r[i])
            atom.radius = elements[symbols[i]]['radius']
            atom.number = elements[symbols[i]]['number']
            atom.id = i % len(self.get_frame_atoms())
            atom.depth = 0
            atom.constrained = False
            if len(ra.constraints) > 0:
                if i in ra.constraints[0].index:
                    atom.constrained = True
            self.queue.append(atom)

    def gfx_queue_box(self):
        try:
            self.get_frame_atoms().cell
        except:
            return
        bx = self.get_frame_atoms().cell
        b = np.array([[0, 0, 0], [bx[1][0], bx[1][1], bx[1][2]], [bx[1][0] +
                     bx[0][0], bx[1][1] + bx[0][1], bx[1][2] + bx[0][2]],
                     [bx[0][0], bx[0][1], bx[0][2]], [bx[2][0], bx[2][1],
                     bx[2][2]], [bx[2][0] + bx[1][0], bx[2][1] + bx[1][1],
                     bx[2][2] + bx[1][2]], [bx[2][0] + bx[1][0] + bx[0][0],
                     bx[2][1] + bx[1][1] + bx[0][1], bx[2][2] + bx[1][2] +
                     bx[0][2]], [bx[2][0] + bx[0][0], bx[2][1] + bx[0][1],
                     bx[2][2] + bx[0][2]]])
        index = [[0, 1], [0, 3], [0, 4], [7, 3], [7, 4], [7, 6], [5, 1], [5, 4],
                 [5, 6], [2, 6], [2, 3], [2, 1]]
        for i in index:
            r1 = b[i[0]]
            r2 = b[i[1]]
            boxsteps = 4
            for l in range(boxsteps):
                self.gfx_queue_line(r1 + (r2 - r1) * float(l) / boxsteps,
                                    r1 + (r2 - r1) * float(l + 1) /
                                    boxsteps, [0, 0, 0])
            
    def gfx_center_atoms(self):
        if self.trajectory is None:
            return
        ra = self.trajectory[0].repeat((int(self.repeatx.get_value()), 
                                        int(self.repeaty.get_value()), 
                                        int(self.repeatz.get_value())))
        r = ra.get_positions()
        minx = min(r[:, 0])                         
        miny = min(r[:, 1])                         
        minz = min(r[:, 2])
        maxx = max(r[:, 0])
        maxy = max(r[:, 1])
        maxz = max(r[:, 2])
        midx = minx + (maxx - minx) / 2
        midy = miny + (maxy - miny) / 2
        midz = minz + (maxz - minz) / 2            
        self.center = np.array([midx, midy, midz])
        self.queue_draw()

    def gfx_draw_axes(self):
        axes = np.identity(3) * 32
        axes[0] = np.dot(self.rotation, axes[0])
        axes[1] = np.dot(self.rotation, axes[1])
        axes[2] = np.dot(self.rotation, axes[2])
        x0 = 48
        y0 = self.height - 48
        self.gfx_draw_line(x0, y0, x0 + axes[0][0], y0 - axes[0][1])
        self.gfx_draw_line(x0, y0, x0 + axes[1][0], y0 - axes[1][1])
        self.gfx_draw_line(x0, y0, x0 + axes[2][0], y0 - axes[2][1])
        font = pango.FontDescription("courier sans 12")
        X = self.area.create_pango_layout("x")
        Y = self.area.create_pango_layout("y")
        Z = self.area.create_pango_layout("z")
        attr = pango.AttrList()
        fg = pango.AttrForeground(65535, 0, 0, 0, -1)
        attr.insert(fg)
        X.set_alignment(pango.ALIGN_CENTER)
        Y.set_alignment(pango.ALIGN_CENTER)
        Z.set_alignment(pango.ALIGN_CENTER)
        X.set_font_description(font)
        Y.set_font_description(font)
        Z.set_font_description(font)
        X.set_attributes(attr)        
        Y.set_attributes(attr)        
        Z.set_attributes(attr)        
        self.pixmap.draw_layout(self.black_gc, int(x0 + axes[0][0]) - X.get_pixel_size()[0] / 2, 
                                int(y0 - axes[0][1]) - X.get_pixel_size()[1], X)
        self.pixmap.draw_layout(self.black_gc, int(x0 + axes[1][0]) - Y.get_pixel_size()[0] / 2, 
                                int(y0 - axes[1][1]) - Y.get_pixel_size()[1], Y)
        self.pixmap.draw_layout(self.black_gc, int(x0 + axes[2][0]) - Z.get_pixel_size()[0] / 2, 
                                int(y0 - axes[2][1]) - Z.get_pixel_size()[1], Z)        
        
            
    def gfx_transform_queue(self):
        for i in range(len(self.queue)):
            q = self.queue[i]
            if q.kind == "atom":
                q.r -= self.center
                q.r = np.dot(self.rotation, q.r)
                q.r += self.translate
                q.depth = q.r[2]
            else:                                   
                q.r1 -= self.center
                q.r2 -= self.center
                q.depth -= self.center
                q.r1 = np.dot(self.rotation, q.r1)
                q.r2 = np.dot(self.rotation, q.r2)
                q.depth = np.dot(self.rotation, q.depth)
                q.r1 += self.translate
                q.r2 += self.translate
                q.depth += self.translate
                q.depth = q.depth[2]


    def gfx_sort_queue(self):
        def cmp_queue(a, b):
            if a.depth > b.depth:
                return 1
            else:
                return -1
        self.queue = sorted(self.queue, cmp_queue)
        

    def gfx_draw_queue(self):
        self.screenatoms = []
        s2 = self.zoombutton.get_value() * 2
        w2 = self.width * 0.5
        h2 = self.height * 0.5
        dr = math.cos(math.pi*0.25)
        for q in self.queue:
            if q.kind == "atom":
                r = q.r
                rad = int(q.radius * self.zoombutton.get_value() * self.radiusbutton.get_value())
                x = int(r[0] * s2 + w2)
                y = int(-r[1] * s2 + h2)
                self.gfx_draw_circle(x, y, rad, q.number)
                if self.frozenbutton.get_active():
                    if q.constrained:
                        self.gfx_draw_line(x-rad*dr, y-rad*dr, x+rad*dr, y+rad*dr)
                        self.gfx_draw_line(x-rad*dr, y+rad*dr, x+rad*dr, y-rad*dr)
                self.screenatoms.append([x, y, rad, q.id])
            else:   
                q.r1[0] = q.r1[0] * s2 + w2
                q.r1[1] = -q.r1[1] * s2 + h2
                q.r2[0] = q.r2[0] * s2 + w2
                q.r2[1] = -q.r2[1] * s2 + h2
                self.gfx_draw_line(q.r1[0], q.r1[1], q.r2[0], q.r2[1])


    def gfx_rot_x(self, theta):
        ct = math.cos(theta)
        st = math.sin(theta)
        m = np.array([[1, 0, 0], [0, ct, -st], [0, st, ct]])
        self.rotation = np.dot(m, self.rotation)
    
    
    def gfx_rot_y(self, theta):
        ct = math.cos(theta)
        st = math.sin(theta)
        m = np.array([[ct, 0, st], [0, 1, 0], [-st, 0, ct]])
        self.rotation = np.dot(m, self.rotation)
    
    
    def gfx_rot_z(self, theta):
        ct = math.cos(theta)
        st = math.sin(theta)
        m = np.array([[ct, -st, 0], [st, ct, 0], [0, 0, 1]])
        self.rotation = np.dot(m, self.rotation)
        

    def gfx_clear(self):
        if self.render_cairo:
            self.cairo_context.set_source_rgb(self.background[0], self.background[1], self.background[2])
            self.cairo_context.rectangle(0, 0, self.width, self.height)
            self.cairo_context.fill()            
        else:
            self.pixmap.draw_rectangle(self.background_gc, True, 0, 0, self.width, self.height)

    def gfx_draw_circle(self, x, y, r, element):
        r = max(1,r)
        if self.render_cairo:
            color = elements[element]['color']
            self.cairo_context.arc(x, y, r, 0, math.pi * 2.0)
            self.cairo_context.set_source_rgb(color[0], color[1], color[2])
            self.cairo_context.fill()
            self.cairo_context.set_source_rgb(0, 0, 0)
            self.cairo_context.arc(x, y, r, 0, math.pi * 2.0)
            self.cairo_context.set_line_width(0.5)
            self.cairo_context.stroke()
        else:
            self.pixmap.draw_arc(self.colors[element], True, x - r, y - r, r * 2, r * 2, 0, 64 * 360)
            self.pixmap.draw_arc(self.black_gc, False, x - r, y - r, r * 2, r * 2, 0, 64 * 360)

    def gfx_draw_line(self, x1, y1, x2, y2):
        if self.render_cairo:
            self.cairo_context.set_source_rgb(0, 0, 0)
            self.cairo_context.move_to(x1, y1)
            self.cairo_context.line_to(x2, y2)
            self.cairo_context.set_line_width(0.5)
            self.cairo_context.stroke()
        else:
            self.pixmap.draw_line(self.black_gc, int(x1), int(y1), int(x2), int(y2))        

#
# DATA ------------------------------------------------------------------------------------------
#

    def data_set(self, data):
        if type(data) != type([]):
            data = [data]
        self.trajectory = data
        self.console.set_item("p", self.trajectory)
        self.moviescale.set_range(0, 1)
        self.moviebox.set_sensitive(False)
        if len(self.trajectory) > 1:
            self.moviebox.set_sensitive(True)
            self.moviescale.set_range(0, len(self.trajectory) - 1)
            self.moviescale.connect("value-changed", lambda w: self.queue_draw())
        self.gfx_center_atoms()
        self.queue_draw()

    def data_read(self, filename):
        data = None
        try:
            f = open(filename, 'r')
            data = []
            while True:
                try:
                    data.append(vasp.read_vasp(f))
                except:
                    f.close()
                    break
            if len(data) < 1:
                raise
        except:
            try:
                data = ase.io.read(filename + "@:")
            except:
                try:
                    data = tsase.io.read_con(filename)
                    if len(data) < 1:
                        raise
                except:
                    try:
                        data = tsase.io.read_bopfox(filename)
                        if len(data) < 1:
                            raise
                    except:
                        print "Failed to load", filename
                        return
        self.data_set(data)
        self.set_title(os.path.abspath(filename))

    def data_write(self, filename):
        if filename.endswith(".con"):
            tsase.io.write_con(filename, self.get_frame_atoms())
        elif filename.endswith(".bx"):
                tsase.io.write_bopfox(filename, self.get_frame_atoms())
        else:
            if '.' not in filename:
                vasp.write_vasp(filename, self.get_frame_atoms())
            else:
                ase.io.write(filename, self.get_frame_atoms())
        self.set_title(os.path.abspath(filename))
        

#
# MAIN ------------------------------------------------------------------------------------------
#

if __name__ == "__main__":
    q = xyz()
    gtk.main()



























