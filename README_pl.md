# Theremino Spectrometer V5.0 – Cross-Platform Python Edition (Linux / macOS / Windows)

Natywny, wieloplatformowy port programu **Theremino Spectrometer V5.0** napisany w Pythonie (**PyQt6 + OpenCV + NumPy**) do obsługi spektrometrów optycznych USB opartych o matrycę CMOS i siatkę dyfrakcyjną (np. popularnych spektrometrów z AliExpress z kontrolerem **Sonix USB 2.0 Camera `0c68:6464` / `USB-ZH`**).

---

## ✨ Główne możliwości

1. **Wieloplatformowa obsługa kamer UVC (Plug & Play):**
   - **Linux:** natywny backend jądra `V4L2` (`cv2.CAP_V4L2`, urządzenia `/dev/video*`).
   - **macOS:** natywny framework `AVFoundation` (`cv2.CAP_AVFOUNDATION`).
   - **Windows:** natywny interfejs `DirectShow` (`cv2.CAP_DSHOW` – identyczny jak w oryginalnym `Theremino_Spectrometer.exe`).
   - **Tryb Demo:** automatyczna symulacja widma świetlówki kompaktowej (CFL Hg/Eu) z szumem matrycy, gdy kamera nie jest podłączona.
2. **Wierny port potoku przetwarzania sygnału (DSP) z Theremino V5.0:**
   - Wybór wycinka obrazu **ROI** w promilach (`Start X`, `End X`, `Start Y`, `Size Y`) oraz odbicie lustrzane `FlipH`.
   - Ekstrakcja luminancji dla każdej kolumny pikseli ($0.299\,R + 0.587\,G + 0.114\,B$) z podglądem paska kamery na żywo i automatycznym wzmocnieniem kontrastu (AGC).
   - **Kalibracja wielopunktowa** długości fali ($\text{nm}$) – pełna kompatybilność z plikami `Calib*.txt` z projektu Theremino oraz wbudowany edytor punktów kalibracyjnych.
   - **Filtry czasowe i przestrzenne:** *Rising Speed*, *Falling Speed*, dwuprzebiegowy filtr wygładzający *Spatial Averaging* oraz skala nieliniowa (*Log Scale*).
   - **Pomiary względne:** odejmowanie tła (*Background / Dark Frame*) oraz dzielenie przez widmo odniesienia (*Reference*) do badań absorpcji i transmisji.
3. **Kolorymetria CIE 1931 i temperatura barwowa (CCT w Kelwinach):**
   - Wyznaczanie współrzędnych chromatycznych **CIE 1931 $(x, y)$** przez całkowanie widma z funkcjami obserwatora standardowego $2^\circ$ ($\bar{x}, \bar{y}, \bar{z}$).
   - Estymacja skorelowanej temperatury barwowej **CCT ($\text{K}$)** (wzory McCamy'ego oraz Hernándeza-Andrésa) oraz odchylenia od ciała doskonale czarnego **$D_{uv}$** (metoda Ohno).
   - Rysowanie na żywo teoretycznej **krzywej promieniowania Plancka** dla wyznaczonej temperatury CCT oraz podgląd barwy źródła (próbka sRGB).
   - Obsługa plików korekcji czułości widmowej matrycy CMOS (`IrradianceCoeffs`, np. `LaoCorrections.txt`).
4. **Fotometria – pomiar natężenia oświetlenia w luksach ($\text{lx}$):**
   - Całkowanie widma z krzywą czułości fotopowej ludzkiego oka **$V(\lambda)$** (maksimum w $555\text{ nm}$) z wizualizacją krzywej $V(\lambda)$ na wykresie.
   - Kompensacja czasu naświetlania ($t_{\text{exp}}$) oraz wzmocnienia ($\text{Gain}$) wraz z ostrzeganiem o przesterowaniu (nasyceniu) matrycy.
   - Jednopunktowa kalibracja wzorcem (np. zewnętrznym luksomierzem) zapisywana do pliku `calib/lux_calib.json`.
5. **Eksport danych:**
   - Zapis widma i parametrów pomiaru (w tym `CCT_Kelvin`, `CIE1931_xy`, `Duv`, `Illuminance_Lux`) do pliku `.csv` w formacie zgodnym z Theremino.

---

## 📁 Struktura projektu

```text
theremino-spectrometer/
├── spectrometer.py                  # Główna aplikacja GUI (PyQt6 + OpenCV + NumPy)
├── README.md                        # Dokumentacja projektu i instrukcja uruchomienia
├── Theremino_Spectrometer_INI.txt   # Oryginalny plik konfiguracyjny Theremino V5.0
├── calib/                           # Pliki kalibracji długości fali (nm) oraz luksomierza
│   ├── Calib1.txt                   # Domyślna kalibracja 4-punktowa (326 nm – 722 nm)
│   ├── Calib2.txt
│   ├── Calib3.txt
│   ├── Calib_WEBCAM.txt
│   └── lux_calib.json               # Współczynnik kalibracji luksomierza (tworzony po kalibracji)
└── irradiance/                      # Pliki korekcji czułości widmowej matrycy CMOS
    ├── LaoCorrections.txt           # Korekcja czułości matrycy z projektu Theremino
    └── Coeffs_FLAT.txt              # Płaska charakterystyka (brak korekcji)
```

---

## 🚀 Instalacja i uruchomienie na każdym systemie

Wymagany jest **Python 3.9+** oraz trzy biblioteki: `opencv-python`, `numpy` oraz `PyQt6`.

---

### 🐧 1. Linux (Ubuntu / Debian / Mint / Fedora / Arch)

#### Krok 1: Instalacja zależności
```bash
sudo apt update
sudo apt install python3 python3-pip python3-opencv python3-numpy python3-pyqt6 acl v4l-utils
```
*(Alternatywnie w wirtualnym środowisku `venv`: `pip install opencv-python numpy PyQt6`).*

#### Krok 2: Uprawnienia do kamery USB (`/dev/video*`)
Podłącz spektrometr do portu USB. Możesz sprawdzić czy system go widzi komendą `lsusb` (powinno pojawić się np. `ID 0c68:6464 Sonix Technology Co., Ltd. USB 2.0 Camera`).

Aby zwykły użytkownik miał dostęp do urządzenia `/dev/video0`:
```bash
# Nadaje uprawnienia natychmiast w bieżącej sesji (bez restartu):
sudo setfacl -m u:$USER:rw /dev/video*

# Dodaje użytkownika na stałe do grupy video (działa po przelogowaniu):
sudo usermod -aG video $USER
```

#### Krok 3: Uruchomienie aplikacji
```bash
python3 spectrometer.py
```

---

### 🍎 2. macOS (Intel oraz Apple Silicon M1 / M2 / M3 / M4)

macOS posiada wbudowany sterownik UVC i wykrywa spektrometr natychmiast po podłączeniu (Plug & Play).

#### Krok 1: Instalacja Pythona i bibliotek
Jeśli używasz menedżera pakietów [Homebrew](https://brew.sh/):
```bash
brew install python
pip3 install opencv-python numpy PyQt6
```

#### Krok 2: Uruchomienie aplikacji
```bash
python3 spectrometer.py
```

#### Uwagi dla macOS:
1. **Zgoda na dostęp do kamery:** Przy pierwszym uruchomieniu macOS wyświetli okno systemowe z pytaniem, czy zezwolić aplikacji *Terminal* (lub *VS Code* / *iTerm*) na dostęp do kamery. Kliknij **Pozwól (Allow)**.
2. **Wybór kamery w MacBookach:** W laptopach z wbudowaną kamerą FaceTime pozycja `Kamera #0` to zazwyczaj kamera nad ekranem. Po podłączeniu spektrometru USB wybierz z rozwijanej listy w lewym górnym rogu programu pozycję **`Kamera #1`** (lub kolejną).

---

### 🪟 3. Windows (Windows 10 / Windows 11)

Windows automatycznie instaluje standardowy sterownik UVC po podłączeniu spektrometru i udostępnia go przez interfejs DirectShow.

#### Krok 1: Instalacja Pythona i bibliotek
1. Pobierz i zainstaluj Python 3 ze strony [python.org](https://www.python.org/downloads/) (podczas instalacji zaznacz opcję **"Add python.exe to PATH"**) lub z Microsoft Store.
2. Otwórz wiersz poleceń (`cmd` lub `PowerShell`) i zainstaluj wymagane pakiety:
   ```powershell
   pip install opencv-python numpy PyQt6
   ```

#### Krok 2: Uruchomienie aplikacji
Przejdź w terminalu do katalogu projektu i wpisz:
```powershell
python spectrometer.py
```

#### Uwagi dla Windowsa:
* Upewnij się, że oryginalny program `Theremino_Spectrometer.exe` (ani żadna inna aplikacja używająca kamery) nie jest uruchomiony w tle, ponieważ w systemie Windows tylko jeden proces naraz może korzystać ze strumienia kamery USB.
* Jeśli w laptopie masz wbudowaną kamerę internetową, przełącz na górnym pasku programu źródło z `Kamera #0` na **`Kamera #1`**.

---

## 🎯 Jak prawidłowo skalibrować spektrometr?

### 1. Ustawienie paska obrazu (ROI)
Na górnym podglądzie z kamery widoczna jest pozioma tęcza rozszczepionego światła oraz pomarańczowa ramka:
* Dopasuj parametry **`Start Y`** oraz **`Size Y`** tak, aby wycinek obejmował najjaśniejszy środkowy pas widma (bez czarnego tła powyżej i poniżej).
* Jeśli kolory tęczy biegną odwrotnie (czerwony po lewej, fioletowy po prawej), zaznacz lub odznacz **`Odwróć X (FlipH)`** – prawidłowo niebieski/fioletowy ($400\text{ nm}$) powinien być po lewej stronie, a czerwony ($700\text{ nm}$) po prawej.

### 2. Kalibracja długości fali ($\text{nm}$) na świetlówce (CFL) lub laserze
1. Skieruj spektrometr na **świetlówkę kompaktową (CFL)** lub ekran oświetlony świetlówką – ma ona bardzo ostre, niezmienne fizycznie piki rtęci i europu:
   - **Niebieski pik rtęci (Hg):** $435.8\text{ nm}$
   - **Zielony pik rtęci (Hg):** $546.1\text{ nm}$
   - **Pomarańczowo-czerwony pik europu (Eu):** $611.6\text{ nm}$
2. Możesz też użyć wskaźnika laserowego (np. fioletowy $405\text{ nm}$, zielony $532\text{ nm}$, czerwony $650\text{ nm}$).
3. Kliknij przycisk **`Edytuj kalibrację...`** i wpisz pozycję piksela odczytaną dla znanych pików oraz odpowiadające im długości fali w $\text{nm}$.

### 3. Pomiar temperatury barwowej CCT ($\text{K}$) oraz luksów ($\text{lx}$)
1. **Odznacz `Auto Ekspozycja`** w panelu *Sterowanie matrycą V4L2*.
2. Ustaw suwak **`Eksp.`** tak, aby wykres badanego źródła światła był wysoki (ok. $60\text{–}85\%$ skali), ale **nie spłaszczony u góry** (brak czerwonego komunikatu o przesterowaniu matrycy).
3. Zasłoń całkowicie wlot spektrometru i kliknij **`Tło (Background)`**, aby odjąć prąd ciemny i szum matrycy.
4. Zaznacz **`Korekcja (LaoCorrections)`**, aby skompensować spadek czułości matrycy CMOS na krańcach pasma widzialnego.
5. Aby wyskalować **Luksomierz ($\text{lx}$)**: połóż obok wlotu spektrometru wzorcowy luksomierz (lub telefon z aplikacją luksomierza), kliknij **`Kalibruj wzorcem (lx)...`** i wpisz odczytaną wartość w luksach. Współczynnik zostanie zapisany automatycznie w `calib/lux_calib.json`.
