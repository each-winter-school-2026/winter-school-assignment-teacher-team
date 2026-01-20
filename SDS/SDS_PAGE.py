
from utils.parseFasta import parseFasta
from _EACH.protein import Protein
from matplotlib import pyplot as plt
import numpy as np
from modules import signal
import base64
from io import BytesIO
from pathlib import Path




# Amino acid molecular weights (in Daltons)
AMINO_ACID_WEIGHTS = {
    #Canonical 20 amino acids
    'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2, 'E': 147.1,
    'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2, 'L': 131.2, 'K': 146.2,
    'M': 149.2, 'F': 165.2, 'P': 115.1, 'S': 105.1, 'T': 119.1, 'W': 204.2,
    'Y': 181.2, 'V': 117.1,
    # Non canonical
    'U': 168.1, 'O': 255.3,
}
# Ambiguous codes (average weights)
AMINO_ACID_WEIGHTS.update({
    'B': (AMINO_ACID_WEIGHTS['D'] + AMINO_ACID_WEIGHTS['N']) / 2,  # Aspartic acid or Asparagine
    'Z': (AMINO_ACID_WEIGHTS['E'] + AMINO_ACID_WEIGHTS['Q']) / 2,  # Glutamic acid or Glutamine
    'J': (AMINO_ACID_WEIGHTS['I'] + AMINO_ACID_WEIGHTS['L']) / 2,  # Isoleucine or Leucine
    'X': sum(AMINO_ACID_WEIGHTS.values()) / len(AMINO_ACID_WEIGHTS)  # Any amino acid
})

AMINO_ACID_HYDROPHOBICITY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'E': -3.5,
    'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
    'Y': -1.3, 'V': 4.2
}


def calculateProteinWeight(sequence):
    """Calculate the molecular weight of a protein sequence.

    Parameters
    ----------
    sequence
        The protein sequence as a string.

    Returns
    -------
        The molecular weight of the protein sequence as a rounded float in kilodaltons.
    """
    return round(sum([AMINO_ACID_WEIGHTS.get(aa,0) for aa in sequence])/1000,ndigits=2)

def convertProteinAbundance(proteins=None,abundances=None):
    #TODO rewrite for clarity
    if proteins is not None:
        abundance = np.array([p.get_abundance() for p in proteins])
    elif abundances is not None:
        abundance = abundances
    else:
        raise ValueError("Either proteins or plottableEvents must be provided.")

    # Minimum and maximum abundance for scaling
    # Minimum opacity is 0 and maximum is 1
    opacityMinimumAbundance = 1e-2
    opacityMaximumAbundance = 1e-0
    minOpacity = 0
    maxOpacity = 1
    
    # Minimum and maximum abundance for line width scaling
    # Minimum line width is ... and maximum is ...
    widthMinimumAbundance = 1e-1
    widthMaximumAbundance = 1e1
    minLineWidth = 1
    maxLineWidth = 25
    
    # Clip abundances to the specified min/max for opacity and width
    ab_op = np.clip(abundance, opacityMinimumAbundance, opacityMaximumAbundance)
    ab_w = np.clip(abundance, widthMinimumAbundance, widthMaximumAbundance)
    print(min(ab_op), max(ab_op), min(ab_w), max(ab_w))

    # Min-max scale to 0-1 with guards against zero division
    den_op = float(opacityMaximumAbundance - opacityMinimumAbundance)
    if den_op == 0:
        scaled_op = np.zeros_like(ab_op, dtype=float)
    else:
        scaled_op = (ab_op - opacityMinimumAbundance) / den_op

    den_w = float(widthMaximumAbundance - widthMinimumAbundance)
    if den_w == 0:
        scaled_w = np.zeros_like(ab_w, dtype=float)
    else:
        scaled_w = (ab_w - widthMinimumAbundance) / den_w

    # Map scaled values to requested opacity and line width ranges
    opacity = minOpacity + (maxOpacity - minOpacity) * scaled_op
    linewidth = minLineWidth + (maxLineWidth - minLineWidth) * scaled_w

    print(min(opacity), max(opacity), min(linewidth), max(linewidth))
    return opacity, linewidth


def binProteinsByWeight(proteins):
    """Bin proteins by their molecular weight.

    Parameters
    ----------
    proteins
        A list of Protein objects to be binned.

    Returns
    -------
        A list of bins, where each bin is a list of proteins.
    """
    if not proteins:
        return []

    weights = np.array([p.get_weight() for p in proteins], dtype=float)
    max_w = float(np.max(weights))

    # If all weights are zero or negative, put everything in one bin
    if max_w <= 0:
        return [list(proteins)]

    bins = np.concatenate((np.arange(0,50,1), np.arange(50,100,2.5),np.arange(100,150,3.5), np.arange(150,261,5)))
    # Assign weights to bins and build result
    idxs = np.digitize(weights, bins) - 1

    p_bins = [[] for _ in range(len(bins))]
    for prot, idx in zip(proteins, idxs):
        p_bins[int(idx)].append(prot)

    return p_bins


def convertBinsToPlottableEvents(bins):
    proteinEvents = []
    for binIndex, binProteins in enumerate(bins):
        proteinWeights = [p.get_weight() for p in binProteins]
        proteinAbundances = [p.get_abundance() for p in binProteins]
        averageWeight = np.average(proteinWeights,weights=proteinAbundances) if proteinWeights else 0
        abundanceSum = sum(proteinAbundances)
        if abundanceSum == 0:
            continue
        proteinEvents.append({
            'binIndex': binIndex,
            'averageWeight': averageWeight,
            'abundanceSum': abundanceSum,
            'proteinCount': len(binProteins),
            'proteins': []#binProteins
        })

    return proteinEvents

# def virtualSDSPage(proteins=Protein.getAllProteins()):
    
#     # Filter proteins with abundance >0
#     proteins = [p for p in proteins if p.get_abundance() > 0]
#     import random
#     random.seed(1)
#     #print(len(proteins))
#     #proteins = random.choices(proteins, k=2000) 
#     binnedProteins = binProteinsByWeight(proteins)
#     plottableEvents = convertBinsToPlottableEvents(binnedProteins)

#     weights = [pe['averageWeight'] for pe in plottableEvents if pe['proteinCount'] > 0]
#     abundances = [pe['abundanceSum'] for pe in plottableEvents if pe['proteinCount'] > 0]
    
#     opacities,widths = convertProteinAbundance(abundances=abundances)
#     fig, ax = plt.subplots(figsize=(10, 2))
#     ax.axes.get_yaxis().set_visible(False)
#     # create a list-of-sequences so each weight is a separate line, and build per-line RGBA colors
#     positions = [[w] for w in weights]
#     print("Plottable events:")
#     for i in plottableEvents:
#         print(i)
#     heights = [len(weights)*1e10 for _ in weights]
#     heights = [1000 for _ in weights]
#     colors = [(0.0, 0.0, 1.0, float(a)) for a in opacities]  # blue with per-line alpha
#     ax.eventplot(positions=positions,linelengths=heights, linewidths=widths, colors=colors, orientation='horizontal')
#     ax.set_title('Virtual SDS-PAGE Protein Weight Distribution')
#     ax.set_xlabel('Molecular Weight (kDa)')
#     ax.set_xscale('log')
#     ax.minorticks_off()  
#     ax.set_xticks([10,15,20,30,40,50,60,80,110,160,260], labels=['10','15','20','30','40','50','60','80','110','160','260'])
#     ax.set_xlim(0,265)
#     #ax.set_ylim(0.4, 1.6)
    
#     plt.tight_layout()
#     #plt.show()
    
    
    
def virtualSDSPage_2DGaussian(proteins=Protein.getAllProteins(),
                                blend='max',  
                                xmin=8, xmax=265,
                                Ny=400, Nx=1200,
                                sigma_x_px=3.0, sigma_y_frac=0.4,
                                amp_range=(0.3, 1.0),
                                ymin_frac=0.3, ymax_frac=0.7,
                                order_y=5,):
    """
        Create and display a schematic "virtual SDS-PAGE" image by rendering each protein
        (or protein bin) as a 2D Gaussian / super-Gaussian intensity blob on an
        image plane that is displayed with matplotlib. The horizontal axis represents
        molecular weight (or migration distance) mapped into image columns with a
        hybrid linear-log mapping; the vertical axis is a normalized lane dimension
        [0, 1] where band shape is controlled by a super-Gaussian vertical profile.
        This function is intended as a visual, schematic representation of
        banding patterns derived from lists of proteins (or bins derived from
        proteins). The function relies on two helper functions defined externally:
        - binProteinsByWeight(proteins) -> returns binned data by molecular weight
        - convertBinsToPlottableEvents(binnedProteins) -> returns a list/dict of
            plottable events; each event must provide at least the keys:
            'averageWeight'  (float): representative molecular weight for the bin
            'abundanceSum'   (float): summed abundance for the bin
            'proteinCount'   (int): number of proteins in the bin (used to filter zeros)
        Important: This function draws the image and shows it with matplotlib but does
        not return the figure or axes objects. If you need to further modify or
        save the figure, call plt.gcf() / plt.gca() after running, or modify the
        function to return (fig, ax).
        Parameters
        ----------
        proteins : sequence
            Iterable of protein-like objects. Each object must implement a method
            .get_abundance() that returns a numeric abundance. Proteins with
            abundance <= 0 are filtered out before plotting. The function expects
            at least one plottable event (non-zero abundance) after binning; otherwise
            numpy operations (e.g., abundances.max()) will raise an exception.
        xmin, xmax : float, float (default: 8, 265)
            The numeric limits of the horizontal axis (representing molecular weight
            or migration scale). These are used both for the final axis extent and to
            shape the mapping from weight to image column. The axis is displayed on a
            log scale via ax.set_xscale('log'), but the mapping from weight to column
            uses a hybrid linear/log mapping controlled by an internal alpha.
        Ny, Nx : int, int (default: 400, 1200)
            Output image dimensions in pixels (rows, columns). Ny controls the
            vertical resolution of the lane (range [0,1]), and Nx controls horizontal
            resolution (mapping of molecular weight to image columns). Higher values
            produce smoother blobs at the cost of memory and CPU.
        sigma_x_px : float (default: 3.0)
            Base horizontal standard deviation (in pixels) used for Gaussian blobs.
            This is scaled per-band according to the band weight (power-law scaling;
            see notes below). Expressed in pixels of the Nx grid.
        sigma_y_frac : float (default: 0.05)
            Vertical spread for the band profile expressed as a fraction of the
            vertical axis height (which is normalized to 1.0). This is the characteristic
            "sigma" for the vertical profile before applying the super-Gaussian exponent.
        amp_range : tuple(float, float) (default: (0.3, 1.0))
            (min_amplitude, max_amplitude) used to linearly map per-bin abundances
            (normalized to the maximum abundance encountered) into the displayed
            amplitude (intensity) of each blob. If all abundances are zero, all blobs
            will be assigned the midpoint of amp_range.
        ymin_frac, ymax_frac : float, float (currently unused; defaults 0.3, 0.7)
            Present in the signature for potential vertical placement control. The
            current implementation uses ycenter = 0.5 (middle of the lane) for all
            bands. These parameters are preserved for future extensions where vertical
            placement varies with abundance or other metadata.
        order_y : int (default: 5)
            The exponent controlling the vertical profile shape. The vertical profile
            uses exp(-0.5 * (|y - 0.5| / sigma_y_frac) ** (2 * order_y)), which produces
            a super-Gaussian when order_y > 1. Higher values of order_y create flatter
            (top-hat-like) bands with sharper edges; order_y = 1 yields a standard
            Gaussian vertical profile.
        blend : {'max', 'sum'} (default: 'max')
            Pixel-wise blending strategy when multiple blobs overlap:
            - 'max' : I(x,y) = max(I_old(x,y), blob(x,y)). This preserves the highest
                        intensity at each pixel (useful to prevent unrealistic
                        saturation where many overlapping bands would sum).
            - 'sum' : I(x,y) = I_old(x,y) + blob(x,y), clipped to the [0,1] range at
                        the end. This enables additive intensity accumulation.
            Other values are treated as 'sum' in the code path (i.e., any non-'max'
            string leads to additive blending).
        Behavior / algorithmic details
        -----------------------------
        1. Filtering:
            - Proteins with abundance <= 0 are removed immediately.
        2. Binning and plottable events:
            - proteins are passed to binProteinsByWeight(), then to
                convertBinsToPlottableEvents(). The latter must yield iterable items (e.g.,
                dicts) with keys 'averageWeight', 'abundanceSum', and 'proteinCount'.
            - Only events with 'proteinCount' > 0 are considered.
        3. Amplitude mapping:
            - Abundances (per-plottable-event) are normalized by the maximum abundance
                observed (abundance / ab_max) to produce values in [0,1]. These are then
                linearly mapped to amp_range.
            - If ab_max == 0 (all zero abundances after binning), all amplitudes are
                set to the midpoint of amp_range.
        4. Horizontal position mapping:
            - A hybrid mapping maps averageWeight -> column index in [0, Nx-1].
                hybrid_map(w) = (1-alpha) * linear_scale + alpha * log_scale with
                alpha = 0.1 (internal constant).
                - linear_scale = (w - xmin) / (xmax - xmin)
                - log_scale   = (log(w) - log(xmin)) / (log(xmax) - log(xmin))
            - The resulting value is clipped to [0,1] and scaled to integer column
                indices across Nx columns. This mapping gives primarily linear placement
                with a slight log bias to compress large-weight spread.
        5. Horizontal band width scaling:
            - The horizontal sigma for each band is computed as:
                sigma_w = sigma_x_px * (w / w_ref) ** beta
                where w_ref = 50.0 and beta = 0.4 (internal constants).
            - Larger molecular weights get proportionally wider horizontal Gaussians,
                according to this power-law scaling.
        6. Band rendering:
            - For each plottable event the code constructs separable 1D profiles:
                gx(x) = exp(-0.5 * ((x - x_center) / sigma_w)**2)
                gy(y) = exp(-0.5 * (| (y - 0.5) / sigma_y_frac | ** (2 * order_y)))
                then constructs an outer product blob = amplitude * gy ⊗ gx.
            - Blobs are either blended by pixelwise maximum or summed into the image.
        7. Color & plotting:
            - The intensity image I (shape Ny x Nx) is converted to an RGB image where
                R = G = 1 - I, B = 1.0. This yields bluish bands on a white-ish background:
                higher I -> more blue; lower I -> near white.
            - A matplotlib figure and axes are created (figsize=(10,1.5)) and the RGB
                image is displayed with origin='lower', aspect='auto', extent=(xmin,xmax,0,1),
                and interpolation='nearest'. The x-axis is set to log scale and a fixed
                set of xticks is applied:
                [10, 15, 20, 30, 40, 50, 60, 80, 110, 160, 260]
                The y-axis is hidden; tight_layout() is applied.
        Return value
        ------------
        None
            The function draws a matplotlib figure/axes and does not explicitly return
            them. If you need programmatic access, call plt.gcf()/plt.gca() after the
            function or modify the function to return (fig, ax).
        Exceptions / edge-cases
        -----------------------
        - If no proteins remain after the initial abundance filter or if convertBinsToPlottableEvents
            yields no events with 'proteinCount' > 0, then the code that computes
            abundances.max() will raise a ValueError. Ensure there is at least one
            non-zero abundance event before calling, or add a guard in caller code.
        - The function expects positive numeric weights; passing non-positive or
            non-finite weights can lead to invalid log operations in the hybrid mapping.
        - Expects numpy as np and matplotlib.pyplot as plt to be imported in the calling
            module namespace and the helper functions binProteinsByWeight and
            convertBinsToPlottableEvents to be defined and available.
        Performance
        -----------
        - Time complexity is O(N_bands * Nx * Ny) in the worst case because for each
            band an Ny-by-Nx outer-product blob is computed (implemented as
            outer(gy, gx) which is Ny*Nx work per band). Memory usage is dominated by
            the Ny x Nx float arrays (I and later RGB).
        - To improve performance for many bands or high resolution:
            - Reduce Ny and/or Nx.
            - Pre-compute gx only for unique sigma_w values or reuse cached blobs for
            repeated/identical parameters.
            - Use lower Ny and a more coarse vertical resolution if vertical detail is
            not critical.
        Usage examples
        --------------
        Basic usage (assuming helper functions are available and proteins is a list):
            virtualSDSPage_2DGaussian(proteins)
        Customize resolution and band shape:
            virtualSDSPage_2DGaussian(proteins, Ny=200, Nx=800, sigma_x_px=4.0, order_y=3)
        Use additive blending (may saturate if many bands overlap):
            virtualSDSPage_2DGaussian(proteins, blend='sum')
        Implementation notes and possible extensions
        -------------------------------------------
        - Currently, vertical centers are fixed at y = 0.5. To introduce vertical
            separation (e.g., to simulate multiple lanes or abundance-dependent vertical
            shifts), replace the fixed ycenter with a per-band ycenter value derived from
            protein metadata or abundance.
        - The hybrid linear-log mapping uses a fixed alpha = 0.1. Expose alpha as a
            parameter for finer control of the angular/compression behavior.
        - Consider returning (fig, ax) for programmatic use, and/or accept an optional
            Axes instance to draw onto.
        - Consider normalizing the final RGB or applying gamma correction for more
            physically plausible contrast.
        Authors / credits
        -----------------
        - Designed as a schematic visualizer for SDS-PAGE style banding. The function
            composes basic image operations with matplotlib for quick visualization of
            abundance vs. molecular-weight distributions.
                            
    """
    if proteins is None:
        img_path = Path('app/EWOKS_Interface/static/img/sdsNoImgSupplied.png')
        with open(img_path, 'rb') as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    
    # Filter proteins with abundance >0
    proteins = [p for p in proteins if p.get_abundance() > 0]
    
    # --- collect weights and abundances from FASTA headers ---
    binnedProteins = binProteinsByWeight(proteins)
    plottableEvents = convertBinsToPlottableEvents(binnedProteins)
    # After line 375 (after computing abundances):
    # Filter plottable events to remove zero-abundance bands
 
    weights = [pe['averageWeight'] for pe in plottableEvents if pe['proteinCount'] > 0]
    abundances = [pe['abundanceSum'] for pe in plottableEvents if pe['proteinCount'] > 0]

    weights = np.array(weights, dtype=float)
    abundances = np.array(abundances, dtype=float)

    # If no plottable bands, produce an empty lane but keep axes
    has_events = weights.size > 0 and abundances.size > 0
    
    

    if has_events:
        # --- map abundances to amplitude range amp_range ---
        # simple linear scaling to [amp_min, amp_max]
        ab_max = abundances.max()
        if ab_max == 0:
            # fallback: all equal amplitude
            amps = np.full_like(abundances, (amp_range[0] + amp_range[1]) / 2.0)
        else:
            ab_norm = abundances / ab_max  # [0,1]
            amps = amp_range[0] + (amp_range[1] - amp_range[0]) * ab_norm
    else:
        # empty image gets filled later; amps unused but keep shape for consistency
        amps = np.array([], dtype=float)

    # === coordinate grids ===
    y = np.linspace(0.0, 1.0, Ny)
    xpix = np.arange(Nx)

    def to_col(w):
        # Simple linear mapping - matplotlib will apply log scaling via set_xscale('log')
        t = (w - xmin) / (xmax - xmin)
        t = np.clip(t, 0, 1)
        return int(np.round(t * (Nx - 1)))

    I = np.zeros((Ny, Nx), dtype=float)

    ycenter = 0.5  # (currently unused; you could vary centers with abundance)

    # === draw each band using abundance-based amplitudes ===
    if has_events:
        for w, A in zip(weights, abundances):
            j = to_col(w)
            xgrid = np.arange(Nx)
            ygrid = y

            w_ref = 50.0
            beta = 0.4
            sigma_w = sigma_x_px * (w / w_ref)**beta

            gx = np.exp(-0.5 * ((xgrid - j) / sigma_w)**2)
            gy = np.exp(-0.5 * (np.abs((y - 0.5) / sigma_y_frac) ** (2 * order_y)))
            blob = A * np.outer(gy, gx)

            if blend == 'max':
                I = np.maximum(I, blob)
            else:
                I += blob

        if blend == 'sum':
            I = np.clip(I, 0, 1)

    RGB = np.empty((Ny, Nx, 3))
    RGB[..., 0] = 1.0 - I
    RGB[..., 1] = 1.0 - I
    RGB[..., 2] = 1.0

    fig, ax = plt.subplots(figsize=(12, 1.5))
    ax.imshow(RGB, origin='lower', aspect='auto',
              extent=(xmin, xmax, 0, 1), interpolation='nearest')
    ax.set_xscale('log')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, 1.0)
    ax.axes.get_yaxis().set_visible(False)
    ax.minorticks_off()
    ax.set_xticks([10,15,20,30,40,50,60,80,110,160,260],
                  labels=['10','15','20','30','40','50','60','80','110','160','260'])
    plt.tight_layout()
    tmpfile = BytesIO()
    plt.savefig(tmpfile, format='png')
    base64EncodedImage = base64.b64encode(tmpfile.getvalue()).decode('utf-8')
    return base64EncodedImage
    

def getExampleSDS_PAGE():
    sequences = parseFasta("data/proteomeHumanUniprot_with_abundance.fasta")

    proteins = []

    for i,hs in enumerate(sequences.items()):
        header, sequence = hs
        proteins.append(Protein(header, sequence))
        
    virtualSDSPage_2DGaussian(proteins, sigma_x_px=3.0, sigma_y_frac=0.4)
    
    tmpfile = BytesIO()
    plt.savefig(tmpfile, format='png')
    
    base64EncodedImage = base64.b64encode(tmpfile.getvalue()).decode('utf-8')
     
    return base64EncodedImage




def main():
    sequences = parseFasta("data/proteomeHumanUniprot_with_abundance.fasta")

    proteins = []

    for i,hs in enumerate(sequences.items()):
        header, sequence = hs
        proteins.append(Protein(header, sequence))
     
       
        
    virtualSDSPage_2DGaussian(proteins)
    Protein.proteinImmunoaffinityDepletion({'ALBU_HUMAN':10.98})
    virtualSDSPage_2DGaussian(proteins)
    Protein.signalPeptideCleavage()
    virtualSDSPage_2DGaussian(proteins)
    plt.show()
    
    #Protein.saveProteinsAsFasta("temp/output")

    #print([calculateProteinWeight(seq) for seq in sequences.values()])


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    print('--- %s seconds ---' % (time.time() - start_time))