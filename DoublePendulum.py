from kivy.app import App                                  # import kelas App dari Kivy — dasar aplikasi Kivy
from kivy.uix.screenmanager import ScreenManager, Screen  # import manajer layar dan kelas layar untuk navigasi antar-screen
from kivy.uix.boxlayout import BoxLayout                  # import layout BoxLayout (susunan hor/ver)
from kivy.uix.gridlayout import GridLayout                # import layout GridLayout (grid)
from kivy.uix.label import Label                          # import widget Label (teks statis)
from kivy.uix.button import Button                        # import widget Button (tombol)
from kivy.uix.slider import Slider                        # import widget Slider (pengatur nilai)
from kivy.uix.textinput import TextInput                  # import widget TextInput (input teks)
from kivy.uix.togglebutton import ToggleButton            # import ToggleButton (tombol toggle/group)
from kivy.graphics import Color, Line, Ellipse, Rectangle, RoundedRectangle, InstructionGroup  # import primitive grafis
from kivy.clock import Clock                              # import Clock untuk scheduling / update berkala
from kivy.core.window import Window                       # import Window untuk konfigurasi jendela (mis. ukuran)
from kivy.metrics import dp                               # import dp (density-independent pixels) untuk ukuran konsisten
from kivy.core.text import Label as CoreLabel             # import CoreLabel untuk menghasilkan tekstur teks custom
from kivy.properties import NumericProperty, BooleanProperty, StringProperty, ListProperty  # import properti untuk binding
from kivy.animation import Animation                      # import Animation untuk efek transisi animasi
import math, time, random                                 # import modul standar: math (matematika), time (waktu), random (acak)

Window.size = (1200, 760)                                 # set ukuran jendela/layar awal aplikasi (width, height)

# ---- Physics integrator (RK4) ----
def derivatives(pendulum, y):                             # fungsi menghitung turunan state untuk sistem double pendulum
    theta1, omega1, theta2, omega2 = y                    # unpack state: sudut dan kecepatan angular untuk dua massa
    m1, m2 = pendulum.m1, pendulum.m2                     # ambil massa dari objek pendulum
    l1, l2 = pendulum.l1, pendulum.l2                     # ambil panjang tali dari objek pendulum
    g = pendulum.g                                        # ambil percepatan gravitasi dari objek pendulum
    delta = theta2 - theta1                               # difference antar sudut (θ2 - θ1)
    cos_d = math.cos(delta); sin_d = math.sin(delta)      # hitung cos dan sin dari delta (dipakai berkali-kali)
    denom1 = (m1 + m2) * l1 - m2 * l1 * cos_d * cos_d     # denominator untuk rumus domega1 (menghindari pembagian nol)
    if abs(denom1) < 1e-12: denom1 = 1e-12 if denom1 >= 0 else -1e-12  # proteksi numeric: kalau hampir nol, ganti nilai kecil
    denom2 = (l2 / l1) * denom1                           # denominator untuk domega2 berdasarkan denom1
    if abs(denom2) < 1e-12: denom2 = 1e-12 if denom2 >= 0 else -1e-12  # proteksi numeric untuk denom2
    domega1 = (                                          # hitung turunan omega1 (akses rumus dinamika)
        m2 * l1 * omega1 * omega1 * sin_d * cos_d +      # bagian inersia/kinetik terkait omega1 dan omega2
        m2 * g * math.sin(theta2) * cos_d +             # bagian gaya gravitasi dari massa kedua yang dikopel
        m2 * l2 * omega2 * omega2 * sin_d -             # kontribusi gerak massa kedua terhadap massa pertama
        (m1 + m2) * g * math.sin(theta1)                # gaya gravitasi pada massa pertama
    ) / denom1                                           # bagi oleh denom1 untuk mendapatkan domega1
    domega2 = (                                          # hitung turunan omega2 (rumus dinamika untuk massa kedua)
        - m2 * l2 * omega2 * omega2 * sin_d * cos_d +   # kontribusi inersia massa kedua
        (m1 + m2) * g * math.sin(theta1) * cos_d -      # gaya gravitasi massa pertama di-coupling ke massa kedua
        (m1 + m2) * l1 * omega1 * omega1 * sin_d -      # kontribusi gerak massa pertama
        (m1 + m2) * g * math.sin(theta2)                # gaya gravitasi langsung pada massa kedua
    ) / denom2                                           # bagi oleh denom2 untuk mendapatkan domega2
    return [omega1, domega1, omega2, domega2]            # kembalikan turunan state sebagai list [dθ1/dt, dω1/dt, dθ2/dt, dω2/dt]

def rk4_step(pendulum, y, dt):                           # fungsi integrator Runge-Kutta 4 untuk satu langkah waktu dt
    k1 = derivatives(pendulum, y)                        # k1 = f(y)
    y2 = [y[i] + 0.5 * dt * k1[i] for i in range(4)]     # sementara y untuk menghitung k2 (y + dt/2 * k1)
    k2 = derivatives(pendulum, y2)                       # k2 = f(y2)
    y3 = [y[i] + 0.5 * dt * k2[i] for i in range(4)]     # sementara y untuk k3
    k3 = derivatives(pendulum, y3)                       # k3 = f(y3)
    y4 = [y[i] + dt * k3[i] for i in range(4)]           # sementara y untuk k4 (y + dt * k3)
    k4 = derivatives(pendulum, y4)                       # k4 = f(y4)
    return [ y[i] + (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) for i in range(4) ]  # rumus kombinasi RK4

# ---- Double pendulum model ----
class DoublePendulum:                                     # kelas model fisika double pendulum: menyimpan state dan update physics
    def __init__(self, m1=1.0, m2=1.0, l1=1.0, l2=1.0, g=9.81):  # konstruktor dengan parameter default
        self.m1, self.m2, self.l1, self.l2, self.g = m1, m2, l1, l2, g  # simpan parameter fisik
        self.state = [math.pi/2, 0.0, math.pi/2, 0.0]    # state awal: θ1, ω1, θ2, ω2 (default berdiri horizontal-ish)
        self.time = 0.0                                  # waktu simulasi (counter)
        self.base_dt = 0.005; self.dt = self.base_dt     # langkah waktu dasar (base_dt) dan dt yang aktif
        self.history = [(0.0, self.state[0], self.state[2])]  # histori waktu & sudut untuk plot/analisis
        self.trail = []                                  # jejak posisi bob kedua untuk digambar (visual trail)
        self.max_history = 5000                          # batas panjang history/trail untuk membatasi memori

    def update(self):                                    # update satu langkah fisika: integrasi + menyimpan jejak
        self.state = rk4_step(self, self.state, self.dt) # integrasikan state menggunakan RK4
        self.time += self.dt                              # tambahkan waktu simulasi
        self.history.append((self.time, self.state[0], self.state[2]))  # simpan waktu dan kedua sudut ke history
        if len(self.history) > self.max_history: self.history.pop(0)    # jaga agar history tidak melebihi max
        theta1, _, theta2, _ = self.state                  # ambil sudut (abaikan omega)
        x1 = self.l1 * math.sin(theta1); y1 = - self.l1 * math.cos(theta1)  # posisi (x1,y1) massa pertama relatif pivot
        x2 = x1 + self.l2 * math.sin(theta2); y2 = y1 - self.l2 * math.cos(theta2)  # posisi (x2,y2) massa kedua
        self.trail.append((x2, y2))                        # tambahkan posisi bob kedua ke trail
        if len(self.trail) > self.max_history: self.trail.pop(0)  # batasi panjang trail

    def get_positions(self):                              # helper untuk mendapatkan posisi pivot, bob1, bob2 (world coords relatif)
        theta1, _, theta2, _ = self.state                 # ambil sudut dari state
        x1 = self.l1 * math.sin(theta1); y1 = - self.l1 * math.cos(theta1)  # hitung posisi bob1
        x2 = x1 + self.l2 * math.sin(theta2); y2 = y1 - self.l2 * math.cos(theta2)  # hitung posisi bob2
        return (0.0, 0.0), (x1, y1), (x2, y2)             # kembalikan pivot, bob1, bob2

# ---- Graph canvas ----
class GraphCanvas(BoxLayout):                             # canvas khusus untuk menggambarkan grafik sudut vs waktu di bagian bawah UI
    def __init__(self, **kwargs):
        super().__init__(**kwargs)                        # panggil konstruktor BoxLayout
        self.size_hint_y = None; self.height = dp(170)    # atur tinggi tetap untuk area grafik
        self.margin_left = dp(60); self.margin_bottom = dp(30)  # margin internal untuk label / ticks
        self.margin_top = dp(12); self.margin_right = dp(16)    # margin atas & kanan
        self.window_duration = 10.0                        # durasi jendela waktu yang ditampilkan (detik)
        self.bg_color = (0.95,0.95,0.95,1)                 # warna latar default (tidak selalu terpakai langsung)
        self.bind(size=self.redraw, pos=self.redraw)      # re-draw saat ukuran/pos berubah

    def redraw(self, *a):                                 # fungsi menggambar ulang seluruh grafik
        self.canvas.clear()                               # bersihkan canvas dulu
        app = App.get_running_app()                       # ambil instance app yang sedang berjalan
        if not hasattr(app,'pendulum') or len(app.pendulum.history) < 2: return  # kalau belum ada data, langsung keluar
        history = app.pendulum.history                    # ambil history dari model pendulum
        t_max = history[-1][0]; t_min = max(0.0, t_max - self.window_duration)  # tentukan rentang waktu yang terlihat
        visible = [(t,a1,a2) for t,a1,a2 in history if t >= t_min]  # filter entri history yang dalam jendela waktu
        if len(visible) < 2:
            visible = history; t_min = visible[0][0]; t_max = visible[-1][0]    # fallback bila sedikit data
        time_range = max(1e-6, t_max - t_min)            # jaga agar time_range tidak nol
        plot_x = self.x + self.margin_left; plot_y = self.y + self.margin_bottom  # koordinat origin plot
        plot_w = max(10, self.width - self.margin_left - self.margin_right)     # lebar area plot
        plot_h = max(60, self.height - self.margin_bottom - self.margin_top)    # tinggi area plot
        all_angles = [a for _,a1,a2 in visible for a in (a1,a2)]  # kumpulkan semua sudut yang terlihat
        max_angle = max(abs(a) for a in all_angles) if all_angles else 1.0      # cari nilai sudut maksimal untuk normalisasi
        if max_angle == 0: max_angle = 1.0                 # proteksi agar tidak membagi nol

        with self.canvas:                                 # mulai blok menggambar
            Color(*App.get_running_app().current_graph_bg)  # set warna background grafik dari tema aplikasi
            Rectangle(pos=self.pos, size=self.size)      # gambar rectangle background untuk seluruh widget
            Color(0.85,0.85,0.85,1)                       # set warna garis grid
            for i in range(1,5):                          # gambar garis horizontal grid
                y_line = plot_y + i*plot_h/5
                Line(points=[plot_x, y_line, plot_x+plot_w, y_line], width=1)
            Color(0.12,0.12,0.12,1)                       # warna untuk axis
            Line(points=[plot_x, plot_y, plot_x, plot_y+plot_h], width=1.5)  # sumbu-y
            Line(points=[plot_x, plot_y, plot_x+plot_w, plot_y], width=1.5)  # sumbu-x

            for i in range(6):                           # gambar tick dan label waktu pada sumbu-x
                t_val = t_min + i * time_range / 5; x_tick = plot_x + i*plot_w/5
                Line(points=[x_tick, plot_y, x_tick, plot_y - dp(6)], width=1.1)  # tick kecil
                lab = CoreLabel(text=f"{t_val:.1f}", font_size=dp(11)); lab.refresh()  # buat tekstur label memakai CoreLabel
                self.canvas.add(Rectangle(texture=lab.texture, pos=(x_tick - lab.texture.size[0]/2, plot_y - dp(26)), size=lab.texture.size))  # gambar label

            for i in range(5):                           # gambar tick & label untuk sumbu-y (angle)
                rel = (i-2)/2.0; angle_val = rel * max_angle
                y_tick = plot_y + i*plot_h/4
                Line(points=[plot_x, y_tick, plot_x - dp(6), y_tick], width=1.1)  # tick kecil horizontal
                lab = CoreLabel(text=f"{angle_val:.2f}", font_size=dp(11)); lab.refresh()  # buat label tekstur
                self.canvas.add(Rectangle(texture=lab.texture, pos=(plot_x - dp(50), y_tick - lab.texture.size[1]/2), size=lab.texture.size))  # gambar label

            labx = CoreLabel(text="Time (s)", font_size=dp(12)); labx.refresh()  # label sumbu x
            self.canvas.add(Rectangle(texture=labx.texture, pos=(plot_x + plot_w/2 - labx.texture.size[0]/2, plot_y - dp(36)), size=labx.texture.size))  # gambar label sumbu x
            laby = CoreLabel(text="Angle (rad)", font_size=dp(12)); laby.refresh()  # label sumbu y
            Push = InstructionGroup()                      # buat instruction group (tidak dipakai penuh di sini, placeholder)
            Push.add(Rectangle(texture=laby.texture, pos=(plot_x - dp(74), plot_y + plot_h/2 - laby.texture.size[1]/2), size=laby.texture.size))  # tambahkan rectangle teks ke group (tidak dipakai rotasi)

            # rotate not strictly necessary; simpler: draw normal left label
            for t,a1,a2 in visible:                        # placeholder loop (tidak melakukan apa-apa di dalamnya)
                pass

            # plot curves
            Color(0.85,0.2,0.2,1)                          # warna garis untuk sudut pertama
            pts1 = []                                      # list titik untuk curve 1
            for t,a1,_ in visible:                         # bangun titik-titik (x,y) untuk curve sudut1
                x = plot_x + ((t - t_min) / time_range) * plot_w
                y = plot_y + plot_h/2 + (a1 / max_angle) * (plot_h/2)
                pts1 += [x,y]
            if len(pts1) >= 4: Line(points=pts1, width=1.6)  # gambar line untuk curve 1 bila cukup titik
            Color(0.15,0.45,0.85,1)                        # warna garis untuk sudut kedua
            pts2=[]                                        # list titik untuk curve 2
            for t,_,a2 in visible:                         # bangun titik-titik untuk curve sudut2
                x = plot_x + ((t - t_min) / time_range) * plot_w
                y = plot_y + plot_h/2 + (a2 / max_angle) * (plot_h/2)
                pts2 += [x,y]
            if len(pts2) >= 4: Line(points=pts2, width=1.6)  # gambar line untuk curve 2 bila cukup titik

    def clear_graph(self):                                # fungsi untuk mereset history grafik
        app = App.get_running_app()                        # ambil instance app
        if hasattr(app,'pendulum'):                        # bila pendulum ada, atur history jadi array awal dengan waktu 0
            app.pendulum.history = [(0.0, app.pendulum.state[0], app.pendulum.state[2])]

# ---- Pendulum visual canvas ----
class PendulumCanvas(BoxLayout):                         # canvas visual utama yang menggambar pendulum (rod + bobs + trail)
    def __init__(self, **kwargs):
        super().__init__(**kwargs)                        # panggil konstruktor BoxLayout
        self.bind(size=self.redraw, pos=self.redraw)      # redraw ketika ukuran/pos berubah
        self.pivot_offset_top = dp(90)                     # offset pivot dari bagian atas widget (visual)
        self.neon_color = (0.12, 0.65, 0.95, 1)            # warna neon untuk trail (RGBA)

    def redraw(self, *a):                                 # fungsi menggambar ulang pendulum
        self.canvas.clear()                               # bersihkan canvas dulu
        app = App.get_running_app()                       # ambil instance app
        if not hasattr(app,'pendulum'): return            # jika belum ada pendulum, keluar
        pend = app.pendulum                                # referensi cepat ke objek pendulum
        _, b1, b2 = pend.get_positions()                  # dapatkan posisi pivot (diabaikan), bob1, bob2
        cx = self.center_x; cy = self.top - self.pivot_offset_top  # hitung pivot world coords (center x dan offset top)
        total_len = max(1e-6, pend.l1 + pend.l2)          # total panjang tali (proteksi nol)
        avail_w = self.width * 0.95; avail_h = (self.height - self.pivot_offset_top) * 0.95  # ruang tersedia untuk menggambar
        avail_r = max(10.0, min(avail_w/2.0, avail_h))    # radius avail untuk skala
        scale = avail_r / total_len                       # skala untuk memetakan panjang fisik ke pixel
        x0,y0 = cx,cy                                     # pivot point di layar
        x1 = cx + b1[0]*scale; y1 = cy + b1[1]*scale      # posisi bob1 dalam koordinat layar (skala applied)
        x2 = cx + b2[0]*scale; y2 = cy + b2[1]*scale      # posisi bob2 dalam koordinat layar

        with self.canvas:                                 # mulai menggambar elemen-elemen
            # background
            Color(*App.get_running_app().current_dark_bg) # warna latar dari tema aplikasi
            Rectangle(pos=self.pos, size=self.size)      # gambar latar widget
            # rods
            Color(0.22,0.22,0.22,1); Line(points=[x0,y0,x1,y1], width=dp(4)); Line(points=[x1,y1,x2,y2], width=dp(4))  # gambar batang pendulum
            # pivot
            Color(0,0,0,1); Ellipse(pos=(x0-dp(5), y0-dp(5)), size=(dp(10),dp(10)))  # gambar pivot sebagai bulatan kecil
            # bobs
            Color(0.9,0.2,0.2,1); Ellipse(pos=(x1-dp(12),y1-dp(12)), size=(dp(24),dp(24)))  # bob pertama (merah-ish)
            Color(0.18,0.45,0.85,1); Ellipse(pos=(x2-dp(12),y2-dp(12)), size=(dp(24),dp(24)))  # bob kedua (biru-ish)
            # trail
            if getattr(App.get_running_app(), 'show_trail', True) and len(pend.trail)>2:  # jika opsi trail aktif dan ada jejak
                Color(self.neon_color[0], self.neon_color[1], self.neon_color[2], 0.18)  # warna trail dengan alpha rendah
                pts=[]                                       # kumpulkan titik trail
                for px,py in pend.trail[-900:]:              # ambil sampai 900 titik terakhir untuk performa
                    pts += [cx + px*scale, cy + py*scale]    # konversi ke koordinat layar
                if len(pts)>=4: Line(points=pts, width=1.3)  # gambar line trail jika cukup titik

# ---- Setup screen (unchanged structure, minimal) ----
class SetupScreen(Screen):                               # layar setup tempat user memilih parameter simulasi
    def __init__(self, **kwargs):
        super().__init__(**kwargs); self.name='setup'     # inisialisasi screen dan beri nama 'setup'
        root = BoxLayout(orientation='horizontal', padding=dp(12), spacing=dp(12))  # root layout horizontal
        form_panel = BoxLayout(orientation='vertical', size_hint=(0.45,1), spacing=dp(8))  # panel kiri untuk form input
        header = Label(text="Simulator - Setup", size_hint=(1,None), height=dp(40))  # header label
        form_panel.add_widget(header)                     # tambahkan header ke panel form
        grid = GridLayout(cols=2, spacing=dp(6), size_hint=(1,None)); grid.bind(minimum_height=grid.setter('height'))  # grid dua kolom untuk label+input
        params = [                                         # daftar parameter (label, default string)
            ("Mass 1 (kg)", "1.00"), ("Mass 2 (kg)", "1.00"),
            ("Length 1 (m)", "1.00"), ("Length 2 (m)", "1.00"),
            ("Gravity (m/s^2)", "9.81"), ("Initial θ1 (rad)", "1.00"),
            ("Initial θ2 (rad)", "2.00")]
        self.inputs={}                                    # dictionary untuk menyimpan pasangan TextInput dan Slider
        for label_text, default in params:                # loop membuat baris input untuk tiap parameter
            grid.add_widget(Label(text=label_text+':', size_hint_y=None, height=dp(32)))  # label parameter
            box = BoxLayout(orientation='horizontal', spacing=dp(6), size_hint_y=None, height=dp(32))  # box untuk slider+textinput
            ti = TextInput(text=default, multiline=False, input_filter='float')  # TextInput untuk memasukkan nilai
            sl_min = -6.28 if "θ" in label_text else 0.0; sl_max = 6.28 if "θ" in label_text else 100.0  # tentukan rentang slider untuk sudut atau non-sudut
            if "Mass" in label_text: sl_min, sl_max = 0.01, 200.0  # jika Mass, ubah rentang slider
            if "Length" in label_text: sl_min, sl_max = 0.01, 10.0  # jika Length, ubah rentang slider
            if "Gravity" in label_text: sl_min, sl_max = 0.01, 50.0  # jika Gravity, ubah rentang slider
            sl = Slider(min=sl_min, max=sl_max, value=float(default), step=0.01)  # buat slider dengan rentang yg ditentukan
            def make_on_slide(ti_ref): return lambda inst,val: setattr(ti_ref,'text',f"{val:.3f}")  # factory: sinkronisasi slider -> TextInput
            sl.bind(value=make_on_slide(ti))                  # bind slider value perubahan ke TextInput
            def make_on_text(sl_ref):                         # factory: sinkronisasi TextInput -> slider
                def fn(inst,val):
                    try:
                        v=float(val); 
                        if v < sl_ref.min: v=sl_ref.min
                        if v > sl_ref.max: v=sl_ref.max
                        sl_ref.value = v
                    except: pass
                return fn
            ti.bind(text=make_on_text(sl))                    # bind TextInput text perubahan ke slider value
            box.add_widget(sl); box.add_widget(ti); grid.add_widget(box)  # tambahkan slider & input ke layout grid
            self.inputs[label_text] = (ti, sl)                # simpan pair di dictionary inputs
        form_panel.add_widget(grid)                          # tambahkan grid ke form_panel
        theme_box = BoxLayout(orientation='horizontal', size_hint=(1,None), height=dp(36), spacing=dp(6))  # box tema
        theme_box.add_widget(Label(text="Theme:", size_hint=(None,1), width=dp(70)))  # label tema
        self.theme_dark = ToggleButton(text="Dark", group="theme", state='down')  # toggle tema Dark default on
        self.theme_light = ToggleButton(text="Light", group="theme")  # toggle Light
        self.theme_blue = ToggleButton(text="Blue", group="theme")    # toggle Blue
        theme_box.add_widget(self.theme_dark); theme_box.add_widget(self.theme_light); theme_box.add_widget(self.theme_blue)  # tambahkan toggles
        form_panel.add_widget(theme_box)                       # tambahkan theme_box ke form_panel
        ctrl_box = BoxLayout(orientation='horizontal', size_hint=(1,None), height=dp(44), spacing=dp(8))  # box kontrol bawah (tombol)
        start_btn = Button(text="Start Simulation", size_hint=(0.55,1)); start_btn.bind(on_press=self.start_simulation)  # tombol Start
        rand_btn = Button(text="Randomize Angles", size_hint=(0.35,1)); rand_btn.bind(on_press=self.randomize_angles)  # tombol random sudut
        ctrl_box.add_widget(start_btn); ctrl_box.add_widget(rand_btn); form_panel.add_widget(ctrl_box)  # tambahkan tombol ke panel
        preview = BoxLayout(orientation='vertical', size_hint=(0.55,1), spacing=dp(8))  # panel preview kanan
        preview.add_widget(Label(text="Preview & Info", size_hint=(1,None), height=dp(36)))  # judul preview
        self.preview_canvas = PendulumCanvas(); preview.add_widget(self.preview_canvas)  # tambahkan canvas preview (visual kecil)
        preview.add_widget(Label(text="Atur parameter di kiri, lalu Start.", size_hint=(1,None), height=dp(40)))  # teks instruksi
        root.add_widget(form_panel); root.add_widget(preview); self.add_widget(root)   # susun root: kiri form, kanan preview, dan tambahkan ke screen

    def randomize_angles(self, inst):                         # fungsi untuk mengacak sudut awal
        a1 = random.uniform(-math.pi, math.pi); a2 = random.uniform(-math.pi, math.pi)  # generate dua sudut acak
        ti1,sl1 = self.inputs["Initial θ1 (rad)"]; ti2,sl2 = self.inputs["Initial θ2 (rad)"]  # ambil textinput & slider terkait
        ti1.text=f"{a1:.3f}"; ti2.text=f"{a2:.3f}"; sl1.value=a1; sl2.value=a2  # set nilai acak ke kedua kontrol

    def start_simulation(self, inst):                         # dipanggil saat user menekan Start Simulation
        app = App.get_running_app()                            # ambil instance aplikasi
        try:
            m1=float(self.inputs["Mass 1 (kg)"][0].text); m2=float(self.inputs["Mass 2 (kg)"][0].text)  # baca m1, m2 dari TextInput
            l1=float(self.inputs["Length 1 (m)"][0].text); l2=float(self.inputs["Length 2 (m)"][0].text)  # baca l1, l2
            g=float(self.inputs["Gravity (m/s^2)"][0].text)  # baca gravitasi
            th1=float(self.inputs["Initial θ1 (rad)"][0].text); th2=float(self.inputs["Initial θ2 (rad)"][0].text)  # baca sudut awal
        except:
            m1,m2,l1,l2,g,th1,th2 = 1.0,1.0,1.0,1.0,9.81,1.0,2.0  # fallback ke nilai default kalau parsing gagal
        app.pendulum = DoublePendulum(m1=m1,m2=m2,l1=l1,l2=l2,g=g)  # buat instance DoublePendulum baru dengan parameter yang dibaca
        app.pendulum.state[0]=th1; app.pendulum.state[2]=th2  # set sudut awal pada state model
        app.pendulum.history=[(0.0,th1,th2)]; app.pendulum.trail=[]  # reset history dan trail
        app.sim_running=True; app.show_trail=True                # set flag simulasi running dan tampilkan trail
        if self.theme_light.state == 'down': app.request_theme('Light')  # jika toggle Light aktif, request theme Light
        elif self.theme_blue.state == 'down': app.request_theme('Blue')  # jika Blue aktif, request theme Blue
        else: app.request_theme('Dark')                         # default: Dark
        app.root.current = 'simulation'                         # pindah layar ke screen 'simulation'

# ---- Simulation screen with drawer ----
class SimulationScreen(Screen):                             # layar simulasi utama dengan drawer kontrol di kanan
    def __init__(self, **kwargs):
        super().__init__(**kwargs); self.name='simulation'   # inisialisasi screen dan beri nama 'simulation'
        root = BoxLayout(orientation='vertical', spacing=dp(8), padding=dp(8))  # root layout vertikal
        # main visual + drawer row
        row = BoxLayout(orientation='horizontal', spacing=dp(8))  # baris horizontal: area visual + drawer
        # left: pendulum visual (big)
        self.left_area = BoxLayout(orientation='vertical', size_hint=(1,1))  # area kiri untuk canvas pendulum
        self.pendulum_canvas = PendulumCanvas()                # instance canvas utama
        self.left_area.add_widget(self.pendulum_canvas)        # tambahkan canvas ke left_area
        row.add_widget(self.left_area)                         # tambahkan left_area ke row
        # drawer (right)
        self.drawer_width_expanded = dp(300); self.drawer_width_collapsed = dp(44)  # ukuran expanded/collapsed drawer
        self.drawer = BoxLayout(orientation='vertical', size_hint=(None,1), width=self.drawer_width_expanded, padding=dp(8), spacing=dp(8))  # container drawer
        # shadow behind drawer (rounded rectangle)
        with self.drawer.canvas.before:
            Color(0,0,0,0.14)                                # warna bayangan
            self._shadow_rect = Rectangle(pos=self.drawer.pos, size=self.drawer.size)  # rectangle yang akan dipindahkan untuk shadow
        self.drawer.bind(pos=self._update_shadow, size=self._update_shadow)  # bind perubahan pos/size drawer untuk update shadow
        # drawer header small (collapse button)
        hdr = BoxLayout(size_hint=(1,None), height=dp(30))    # header kecil untuk drawer
        self.toggle_btn = Button(text="≡", size_hint=(None,1), width=dp(40))  # tombol toggle untuk collapse/expand
        self.toggle_btn.bind(on_press=self.toggle_drawer); self._add_fade(self.toggle_btn)  # bind aksi toggle + tambahkan efek fade
        hdr.add_widget(Label()); hdr.add_widget(self.toggle_btn)  # tambahkan toggle ke header (dengan spacer label)
        self.drawer.add_widget(hdr)                            # tambahkan header ke drawer
        # buttons grid
        grid = GridLayout(cols=2, spacing=dp(8), size_hint=(1,None), height=dp(220))  # grid tombol di drawer
        self.side_resume = Button(text="RESUME"); self._add_fade(self.side_resume); self.side_resume.bind(on_press=self.toggle_pause)  # resume/pause
        self.side_restart = Button(text="RESTART"); self._add_fade(self.side_restart); self.side_restart.bind(on_press=self.reset_sim)  # restart
        self.side_back = Button(text="BACK TO SETUP"); self._add_fade(self.side_back); self.side_back.bind(on_press=self.back_to_setup)  # kembali ke setup
        self.side_trail = Button(text="TRAIL: OFF"); self._add_fade(self.side_trail); self.side_trail.bind(on_press=self.toggle_trail)  # toggle trail on/off
        self.side_clear = Button(text="CLEAR GRAPH"); self._add_fade(self.side_clear); self.side_clear.bind(on_press=self.clear_graph)  # clear graph
        self.side_theme = Button(text="THEME"); self._add_fade(self.side_theme); self.side_theme.bind(on_press=self.cycle_theme)  # cycle theme button
        grid.add_widget(self.side_resume); grid.add_widget(self.side_restart)  # tambahkan tombol ke grid (row1)
        grid.add_widget(self.side_back); grid.add_widget(self.side_trail)     # row2
        grid.add_widget(self.side_clear); grid.add_widget(self.side_theme)    # row3
        self.drawer.add_widget(grid)                           # tambahkan grid tombol ke drawer
        # time display
        timebox = BoxLayout(orientation='vertical', size_hint=(1,None), height=dp(100))  # box untuk menampilkan waktu simulasi
        timebox.add_widget(Label(text="TIME", size_hint=(1,None), height=dp(18)))  # label TIME
        self.big_time = Label(text="00:00:00", font_size=dp(20)); timebox.add_widget(self.big_time)  # label besar untuk menampilkan waktu
        self.drawer.add_widget(timebox)                        # tambahkan timebox ke drawer
        # speed slider + value (value below)
        speedbox = BoxLayout(orientation='vertical', size_hint=(1,None), height=dp(120))  # box untuk speed control
        speedbox.add_widget(Label(text="SPEED", size_hint=(1,None), height=dp(20)))  # label SPEED
        self.speed_slider = Slider(min=0.25, max=10.0, value=1.0, step=0.05)  # slider kecepatan simulasi
        self.speed_slider.bind(value=self.on_speed_change); self._add_fade(self.speed_slider)  # bind perubahan nilai ke handler
        speedbox.add_widget(self.speed_slider)                # tambahkan slider ke speedbox
        self.speed_value = Label(text="1.00x", size_hint=(1,None), height=dp(20)); speedbox.add_widget(self.speed_value)  # label menunjukkan nilai speed
        self.drawer.add_widget(speedbox)                       # tambahkan speedbox ke drawer
        row.add_widget(self.drawer)                            # tambahkan drawer ke row (kanan)
        root.add_widget(row)                                   # tambahkan row (visual + drawer) ke root vertikal
        # graph at bottom
        self.graph = GraphCanvas()                             # instance GraphCanvas di bagian bawah
        root.add_widget(self.graph)                            # tambahkan graph ke root (di bawah)
        self.add_widget(root)                                  # tambahkan root layout ke screen

        # state
        self.paused = False; self.control_open = True          # inisialisasi state paused dan apakah control terbuka
        self.update_ev = Clock.schedule_interval(self.update, 1.0/60.0)  # schedule update dipanggil 60x per detik

    def _update_shadow(self, instance, value):                # callback untuk memperbarui posisi/size shadow rectangle
        try:
            self._shadow_rect.pos = (self.drawer.x - dp(8), self.drawer.y)  # set posisi shadow mengikuti drawer
            self._shadow_rect.size = (self.drawer.width + dp(8), self.drawer.height)  # set ukuran shadow sesuai drawer
        except Exception:
            pass                                              # aman-aman saja bila ada error kecil (ignore)

    def _add_fade(self, widget):                              # helper menambahkan efek fade singkat saat tombol ditekan
        def on_press(inst):
            anim = Animation(opacity=0.55, duration=0.06) + Animation(opacity=1.0, duration=0.12)  # animasi fade-out then in
            anim.start(inst)                                # jalankan animasi pada widget yang ditekan
        widget.bind(on_press=on_press)                       # bind event on_press ke fungsi on_press

    def toggle_drawer(self, instance):                        # fungsi untuk toggle open/collapse drawer
        target_w = self.drawer_width_collapsed if self.control_open else self.drawer_width_expanded  # target lebar tergantung state
        anim = Animation(width=target_w, d=0.22, t='out_cubic')  # buat animasi perubahan lebar
        anim.start(self.drawer)                              # jalankan animasi pada drawer
        self.control_open = not self.control_open            # invert state control_open
        # change button arrow
        self.toggle_btn.text = "›" if self.control_open else "≡"  # ubah teks tombol untuk merefleksikan state

    def toggle_pause(self, instance):                         # fungsi toggle pause/resume simulasi
        self.paused = not self.paused                         # invert paused
        # update labels
        self.side_resume.text = "RESUME" if self.paused else "PAUSE"  # atur teks tombol sesuai state (note: inverted logic naming)

    def toggle_trail(self, instance):                         # fungsi toggle menampilkan trail jejak bob kedua
        app = App.get_running_app()
        app.show_trail = not getattr(app, 'show_trail', True)  # invert flag show_trail pada app
        self.side_trail.text = "TRAIL: ON" if app.show_trail else "TRAIL: OFF"  # update teks tombol

    def on_speed_change(self, instance, value):                # handler saat slider speed berubah
        App.get_running_app().speed = float(value)             # set property speed di app (dipakai untuk mengatur dt)
        self.speed_value.text = f"{value:.2f}x"                # update label yang menampilkan multiplier speed

    def cycle_theme(self, instance):                           # fungsi untuk berganti theme dalam urutan
        app = App.get_running_app()
        if app.theme_name == 'Dark': app.request_theme('Light')  # Dark -> Light
        elif app.theme_name == 'Light': app.request_theme('Blue')  # Light -> Blue
        else: app.request_theme('Dark')                         # Blue -> Dark

    def update(self, dt):                                      # dipanggil tiap frame oleh Clock.schedule_interval
        app = App.get_running_app()
        if hasattr(app,'pendulum'):
            app.pendulum.dt = app.pendulum.base_dt * getattr(app,'speed',1.0)  # atur dt model berdasarkan multiplier speed
        if getattr(app,'sim_running',False) and not self.paused and hasattr(app,'pendulum'):
            sub = max(1, int(max(1.0, getattr(app,'speed',1.0)) * 2))  # tentukan berapa sub-step update untuk kestabilan saat speed tinggi
            for _ in range(sub): app.pendulum.update()              # jalankan beberapa update fisika per frame bila perlu
        if hasattr(app,'pendulum'):
            self.big_time.text = time.strftime("%H:%M:%S", time.gmtime(app.pendulum.time))  # tampilkan waktu simulasi dalam format HH:MM:SS (GMT based)
        # update canvases' background via app animated colors
        try:
            self.pendulum_canvas.redraw(); self.graph.redraw()  # panggil redraw untuk canvas visual dan graph
        except Exception:
            pass                                              # jika terjadi error saat redraw, ignore agar UI tidak crash

    def reset_sim(self, instance):                             # fungsi restart/reset state dinamika (tetap pada sudut awal)
        app = App.get_running_app()
        if hasattr(app,'pendulum'):
            p = app.pendulum; p.state[1]=0.0; p.state[3]=0.0; p.time=0.0  # set omega ke 0 dan waktu ke 0
            p.history=[(0.0,p.state[0],p.state[2])]; p.trail=[]           # reset history dan trail
            self.graph.clear_graph()                        # clear graph juga

    def back_to_setup(self, instance):                         # kembali ke layar setup
        App.get_running_app().root.current = 'setup'           # ubah screen manager ke 'setup'

    def clear_graph(self, instance):                           # clear grafik (panggil fungsi GraphCanvas.clear_graph)
        self.graph.clear_graph()

# ---- Main App ----
class DoublePendulumApp(App):                               # kelas utama aplikasi Kivy
    speed = NumericProperty(1.0)                             # property speed (float) dengan binding Kivy
    show_trail = BooleanProperty(True)                       # property boolean untuk menampilkan trail
    theme_name = StringProperty("Dark")                      # nama theme saat ini
    pendulum = None                                          # placeholder untuk instance DoublePendulum
    sim_running = False                                      # flag apakah simulasi sedang berjalan

    current_dark_bg = ListProperty([0.06,0.06,0.06,1])       # properti warna latar gelap saat ini (rgba)
    current_graph_bg = ListProperty([0.95,0.95,0.95,1])      # properti warna background graph saat ini
    current_light_bg = ListProperty([1,1,1,1])               # properti warna latar terang saat ini

    THEMES = {                                               # definisi palet tema yang tersedia
        'Dark': {'dark_bg':(0.06,0.06,0.06,1), 'light_bg':(1,1,1,1), 'graph_bg':(0.95,0.95,0.95,1)},
        'Light':{'dark_bg':(1,1,1,1), 'light_bg':(1,1,1,1), 'graph_bg':(0.97,0.97,0.97,1)},
        'Blue': {'dark_bg':(0.02,0.03,0.08,1), 'light_bg':(0.02,0.03,0.08,1), 'graph_bg':(0.03,0.05,0.1,1)}
    }

    def build(self):                                         # dipanggil ketika aplikasi dibangun (sebelum on_start)
        self.title = "Double Pendulum - Drawer UI"          # set judul window aplikasi
        sm = ScreenManager(); self.setup = SetupScreen(); self.sim = SimulationScreen()  # buat screen manager dan screen
        sm.add_widget(self.setup); sm.add_widget(self.sim)   # tambahkan screen ke screen manager
        pal = self.THEMES.get(self.theme_name, self.THEMES['Dark'])  # ambil palet warna awal berdasarkan theme_name
        self.current_dark_bg = list(pal['dark_bg']); self.current_graph_bg = list(pal['graph_bg'])  # set warna awal
        self.current_light_bg = list(pal['light_bg'])
        return sm                                           # kembalikan root widget (screen manager) sebagai UI aplikasi

    def on_start(self):                                     # dipanggil saat aplikasi selesai build dan mulai berjalan
        self.pendulum = DoublePendulum(); self.sim_running = False; self.speed = 1.0; self.theme_name='Dark'  # inisialisasi model dan flags

    def request_theme(self, new_theme, duration=0.25):      # fungsi untuk mengubah tema dengan transisi animasi
        if new_theme not in self.THEMES: return             # jika tema tidak valid, abaikan
        start_dark = list(self.current_dark_bg); end_dark = list(self.THEMES[new_theme]['dark_bg'])  # tetapkan warna awal & target untuk dark_bg
        start_graph = list(self.current_graph_bg); end_graph = list(self.THEMES[new_theme]['graph_bg'])  # untuk graph_bg
        start_light = list(self.current_light_bg); end_light = list(self.THEMES[new_theme]['light_bg'])  # untuk light_bg
        steps = max(6, int(60 * duration)); step = {'i':0}  # hitung jumlah langkah interpolasi berdasarkan durasi animasi
        def stepper(dt):
            i = step['i'] + 1; t = i / steps               # progres t di antara 0..1
            def lerp(a,b,t): return a + (b-a) * t          # linear interpolation helper
            self.current_dark_bg = [lerp(start_dark[j], end_dark[j], t) for j in range(4)]  # interpolasi warna dark_bg
            self.current_graph_bg = [lerp(start_graph[j], end_graph[j], t) for j in range(4)]  # interpolasi warna graph_bg
            self.current_light_bg = [lerp(start_light[j], end_light[j], t) for j in range(4)]  # interpolasi warna light_bg
            step['i'] = i
            if i >= steps:
                self.current_dark_bg = end_dark; self.current_graph_bg = end_graph; self.current_light_bg = end_light  # pastikan set akhir persis target
                self.theme_name = new_theme                  # set nama theme ke yang baru
                return False                                # hentikan schedule
            return True                                     # lanjutkan schedule
        Clock.schedule_interval(stepper, duration/steps)      # schedule stepper untuk dijalankan berkala hingga selesai

if __name__ == '__main__':                                 # standar Python: jika file ini dijalankan langsung
    DoublePendulumApp().run()                              # buat instance App dan jalankan aplikasi