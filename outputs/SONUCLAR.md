# Pourazad-style per-unit VSM-GFM contraction model

Aynı MATLAB dosyası güncellendi. H=3.4, K=20, Kd=5, D=25; E=U1=1 ve Xeq=0.3. Ana nominal hız omega0_pu=1. Fiziksel taban 50 Hz ve omega_base=314.159265 rad/s yalnızca referans bilgisidir.

## Model ve normalizasyon

x=[gamma;nu]; gamma_dot=nu-1; nu_dot=[(Pstar-Pe)/nu-D(nu-1)]/(2H); Pe=(E U1/Xeq)sin(gamma). 1/nu korunur. Her gamma0 için nu0=1 ve Pstar=Pe(gamma0).

Bu, istenen makale ölçeğindeki analitik modeldir. Fiziksel EMT modelinde gamma_dot=omega_base(nu-1) olur. Yalnızca açı denkleminden omega_base kaldırılırken frekans denklemi aynı tutulduğu için bu değişiklik fiziksel modelin salt koordinat veya ortak zaman dönüşümü değildir. Grafik zamanı fiziksel saniye olarak yorumlanmaz.

A21=-Pe_gamma/(2H nu); A22=[-(Pstar-Pe)/nu^2-D]/(2H); B2=-E sin(gamma)/(2H Xeq nu). Dengede alpha2=-25/6.8=-3.676470588235. beta_dot_paper=-E cos(gamma0)/(2H Xeq); beta_dot_paper/beta=cot(gamma0).

Ana metrik bu makale katsayısıyla her denge için ayrı tasarlanır ve referans boyunca sabittir. Tam doğrudan F hesabıyla Eq15c/2 uyuşur. Gerçek durum türevi beta_gamma*(nu-1)+beta_nu*nu_dot ayrıca hesaplanır; dengede sıfırdır ve onun farklı metriği sertifika vermez.

Candidate contraction-oriented retuning; optimallik iddiası yoktur. D/(2H)=3.67647058824; D/H=7.35294117647.

## Yeni sonuçlar

- 15.3225991806° < gamma0 < 21.0916682548°; sınır güçleri 0.88084493, 1.19953713 pu.
- 48.6263618306° < gamma0 < 81.4159628061°; sınır güçleri 2.50138420, 3.29599335 pu.

Her sınırın sebebi boundary_causes.csv dosyasında; negatif köşegenler tek başına kullanılmadı, determinant ve tam özdeğer koşulu kontrol edildi.

| Case | gamma0 (deg) | Pstar (pu) | lambda max paper |
|---|---:|---:|---:|
| Below R1 - 12 deg | 12.000000000 | 0.693038969 | 1.749865311 |
| Inside R1 - nominal 1 pu | 17.457603000 | 0.999999993 | -0.131012073 |
| Between R1/R2 - 35 deg | 35.000000000 | 1.911921455 | 0.110107334 |
| Inside R2 - 65 deg | 65.000000000 | 3.021025957 | -0.051357271 |
| Above R2 - 120 deg | 120.000000000 | 2.886751346 | 0.652064080 |

Nominal 1 pu güç açısı 17.457603124 derecedir. Güçler ve 1 pu sınırı ayrıca grafiklenmiştir.

Nominal 1 pu sertifikası: 1; yerel kararlılık: 1; lambda_max(SymF)=-0.131012071369.

Önceden doğrulanmış eski PU ayarı: H=0.1, K=20, Kd=40, D=60; aralık [89.673050904,89.874923488] derece; nominal sertifikasız. Bu bilgi yalnızca karşılaştırmadır. Toplam aralık genişliği değişimi +38.356797465 derece; nominal lambda_max değişimi -337.248545933.

## Doğrulama ve yörüngeler

Ana model Jacobianı, giriş Jacobianı, tam beta türevleri, PU mutlak/sapma koordinat eşdeğerliği ve tüm denge kalıntıları test edildi. Sertifika eşdeğerliği (s11<0, s22<0, det>0) ile lambda_max<0 arasında tüm taramada kontrol edildi.

12 farklı denge noktası matrislerden sınıflandırıldı; tüm matrisler MAT dosyalarında, sınıflandırma 12_point_gfm_classification.csv dosyasında saklanır. Beş temsilci gamma(0)=gamma0+1 derece ve nu(0)=1 ile başlatılır. Bu ortak test bozuntusu bize ait görselleştirme tercihidir; referans makaleden alınan bir büyüklük değildir. Figure-5 tarzı bozuntu şiddeti analizi daha sonra yapılacaktır. Ana figürde 120 derece için özgün dengeden ilk 20 derece ayrılmaya kadar aynı zaman indeksleriyle faz düzlemi ve L(t), ek 10_above_R2_full_time_diagnostic figüründe tam simülasyon gösterilir. L(t) özgün dengeye göredir.

- Below R1 - 12 deg: Returns to original equilibrium; final gamma=12.000000318 derece; güç kalıntısı=+1.80921e-08 pu.
- Inside R1 - nominal 1 pu: Returns to original equilibrium; final gamma=17.457603318 derece; güç kalıntısı=+1.76437e-08 pu.
- Between R1/R2 - 35 deg: Returns to original equilibrium; final gamma=35.000000317 derece; güç kalıntısı=+1.51223e-08 pu.
- Inside R2 - 65 deg: Returns to original equilibrium; final gamma=65.000000317 derece; güç kalıntısı=+7.7848e-09 pu.
- Above R2 - 120 deg: Approaches another equilibrium branch; final gamma=419.999050275 derece; güç kalıntısı=-2.76268e-05 pu.
Eski SI sonuçları outputs/diagnostics altında OLD SI-coordinate / mixed-normalization diagnostic; NOT used as thesis baseline result etiketiyle korunmuştur. Ana çıktılar yalnızca yeni PU hesabına aittir.
