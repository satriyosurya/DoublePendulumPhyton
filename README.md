# DoublePendulumPhyton

## ✨ Fitur Utama

* **Simulasi Fisika Akurat:** Menggunakan integrasi numerik **Runge-Kutta Orde 4 (RK4)** untuk menghitung dinamika *double pendulum* yang terkenal kacau (*chaotic*).
* **Antarmuka Kivy Interaktif:** Antarmuka pengguna (UI) yang responsif dengan panel kontrol *drawer* yang dapat dilipat (*collapsible*).
* **Pengaturan Parameter Dinamis:** Pengguna dapat mengatur Massa, Panjang Tali, Gravitasi, dan Sudut Awal melalui *Slider* dan *Text Input*.
* **Visualisasi Waktu Nyata:** Menampilkan pergerakan pendulum secara visual, termasuk jejak (*trail*) bob kedua.
* **Plot Sudut vs Waktu:** Grafis interaktif di bagian bawah untuk memvisualisasikan $\theta_1$ dan $\theta_2$ terhadap waktu simulasi.
* **Kontrol Simulasi:** Fitur Jeda/Lanjutkan (*Pause/Resume*), Atur Ulang (*Restart*), dan pengatur Kecepatan Simulasi.
* **Tema:** Dukungan untuk berganti tema (Dark, Light, Blue) dengan transisi warna yang mulus.

## ⚙️ Instalasi dan Menjalankan Proyek

Proyek ini membutuhkan Python 3 dan *framework* Kivy.

### Persyaratan
Pastikan Anda telah menginstal Python 3 dan *pip*.

### Langkah-Langkah Instalasi

1.  **Clone Repository:**
    ```bash
    git clone [https://github.com/UsernameAnda/DoublePendulumKivy.git](https://github.com/UsernameAnda/DoublePendulumKivy.git)
    cd DoublePendulumKivy
    ```

2.  **Instal Dependensi:**
    ```bash
    pip install kivy
    ```

3.  **Jalankan Aplikasi:**
    ```bash


## 🛠️ Struktur Kode

Seluruh simulasi dan UI diimplementasikan dalam satu file `DoublePendulum.py`, yang dibagi menjadi beberapa bagian logis:

1.  **Physics Integrator (RK4):** Fungsi `derivatives` dan `rk4_step` yang mengimplementasikan integrasi numerik untuk persamaan diferensial *double pendulum*.
2.  **DoublePendulum Model:** Kelas yang menyimpan state fisika ($\theta_1, \omega_1, \theta_2, \omega_2$) dan mengelola *history* serta *trail*.
3.  **GraphCanvas:** Widget Kivy untuk menggambar plot sudut terhadap waktu.
4.  **PendulumCanvas:** Widget Kivy untuk visualisasi batang, bob, dan jejak (*trail*) pendulum.
5.  **SetupScreen:** Layar untuk mengatur parameter awal simulasi.
6.  **SimulationScreen:** Layar utama simulasi dengan kanvas visual, grafik, dan panel kontrol *drawer*.
7.  **DoublePendulumApp:** Kelas utama aplikasi Kivy, mengelola *ScreenManager*, *Themes*, dan *global properties* (seperti *speed*).
    python DoublePendulum.py
    ```
    *Catatan: Aplikasi dirancang untuk ukuran jendela awal 1200x760.*
