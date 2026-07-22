# Layered SVG authoring guide

This guide defines the SVG contract used by the human multiregion AAV map.
The current application uses a hybrid implementation: named Reactome organ
groups are reused as alpha masks for exact organ silhouettes, while model-only
subregions (brain lobes, kidney medulla, muscle territories, marrow, heart, and
spleen) are supplied by aligned local vector layers. The long-term target is one
authored SVG in which every model region is a named, closed vector path.

## 0. Current runtime mapping

`app/organ-heatmap.tsx` maps the following Reactome fragments directly to ODE
regions. The original gradients remain visible as a low-opacity reference;
the mask receives the live ODE color above it.

| Reactome fragment | Runtime ODE region |
| --- | --- |
| `R-ICO-013680` | six clipped brain subregions |
| `R-ICO-013935` | left and right lung |
| `R-ICO-012959` | liver |
| `R-ICO-012931` | left/right kidney cortex plus medulla inset |
| `R-ICO-013406`, `R-ICO-012904`, `R-ICO-013151` | combined gut |
| `BG` | body silhouette for systemic and muscle territories |

The colored geometry and transparent hit targets are separate. This prevents a
large masked rectangle from intercepting clicks outside the visible organ and
allows keyboard focus for all 24 model regions.

## 1. Document setup

- Canvas/viewBox: `0 0 600 1000`.
- Keep the existing frontal human figure aligned to this coordinate system.
- Work in a single SVG document.
- Use four top-level layers in this order:
  1. `reference`: locked anatomical reference, not interactive.
  2. `regions`: closed shapes that React colors.
  3. `vessels`: injection and circulation paths.
  4. `labels`: optional editing labels; hide or remove before delivery.
- Draw with neutral fills while editing. React replaces region fills at runtime.
- Avoid bitmap effects, blur, and text that depends on a locally installed
  font. A final self-contained authored SVG should not require external masks;
  the current Reactome-backed runtime intentionally uses same-origin fragment
  masks as an intermediate integration strategy.

Recommended structure:

```xml
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 600 1000">
  <g id="reference" pointer-events="none">...</g>
  <g id="regions">
    <path id="brain-frontal" data-region="brain_frontal" d="..." />
    <g id="kidney-left-cortex" data-region="kidney_left_cortex">
      <path d="..." />
      <path d="..." />
    </g>
  </g>
  <g id="vessels" fill="none">...</g>
</svg>
```

Use `data-region` for the model key. The visual `id` may use hyphens, but every
`data-region` value must exactly match the underscore names below.

## 2. Required region IDs

### CNS

- `brain_frontal`
- `brain_parietal`
- `brain_temporal`
- `brain_occipital`
- `brain_deep_gray`
- `brain_cerebellum`
- `brainstem_spinal`

### Thoracic and abdominal organs

- `lung_left`
- `lung_right`
- `heart`
- `liver`
- `spleen`
- `gut`

### Kidney subregions

- `kidney_left_cortex`
- `kidney_left_medulla`
- `kidney_right_cortex`
- `kidney_right_medulla`

The medulla paths should sit above the cortex paths. Do not cut the medulla out
of the cortex unless you intentionally want the cortex to exclude it visually.

### Muscle and systemic tissues

- `muscle_injected_arm`
- `muscle_contralateral_arm`
- `muscle_trunk`
- `muscle_legs`
- `skin_adipose`
- `bone_marrow`
- `rest`

The patient's anatomical left arm appears on the viewer's right in a frontal
view. That side must be `muscle_injected_arm` in the current left-arm IV demo.

## 3. Illustrator workflow

1. Open the current anatomical SVG and set the artboard to 600 by 1000.
2. In the Layers panel create `reference`, `regions`, `vessels`, and `labels`.
3. Move the original illustration into `reference`, reduce opacity to about
   25%, and lock the layer.
4. In `regions`, use the Pen or Curvature tool to trace each required region.
5. Close every shape. Open paths cannot receive a reliable interactive fill.
6. Use Shape Builder or Pathfinder Unite when one model region consists of
   several visual pieces. Alternatively, keep the pieces inside one named group.
7. For nested anatomy, place the outer region first and inner region later. For
   example, kidney cortex below kidney medulla and brain lobes below deep gray.
8. In the Layers panel rename each object or group using the visual hyphenated
   ID, for example `brain-frontal`.
9. Add the exact `data-region` attributes after export in a text editor, or use
   an SVG metadata-capable plugin. Illustrator object names normally become IDs
   but do not automatically create `data-region`.
10. Draw blood vessels as stroked open paths in `vessels`; these do not need to
    be closed. Keep the injection point as a separate circle.
11. Export with SVG 1.1, Styling set to Presentation Attributes, Fonts converted
    to outlines or labels removed, Images embedded, Object IDs set to Layer
    Names, and 3-4 decimal places.
12. Disable responsive scaling only if Illustrator tries to remove the viewBox.
    The delivered SVG must retain `viewBox="0 0 600 1000"`.

## 4. Inkscape workflow

1. Open Document Properties and set the page to 600 by 1000 px.
2. Create the four layers with Layer > Add Layer.
3. Lock the `reference` layer and trace regions with the Bezier tool.
4. Use Path > Object to Path for ellipses or rectangles that should become
   final paths.
5. Use Path > Union for pieces belonging to one region and Path > Difference
   only when a real anatomical hole is needed.
6. Open Object Properties for each path/group and set its ID.
7. Save as Plain SVG. Inkscape SVG contains extra editor metadata that is not
   needed by React.
8. Add or verify `data-region` attributes in the XML Editor or a text editor.

## 5. Boundaries and overlap

- Adjacent areas may share an edge, but avoid visible gaps larger than 1 px.
- Avoid two large opaque regions covering each other. Nested anatomy such as
  medulla inside cortex is the exception.
- Use groups for disconnected pieces that share one ODE state, such as left and
  right long-bone marrow, rather than duplicating a `data-region` value on many
  unrelated top-level objects.
- A frontal body cannot show every deep structure clearly. Use a separate CNS
  inset or kidney inset while retaining the same region IDs.
- Do not encode scientific values as permanent SVG colors. Color is supplied by
  the React scale from ODE output.

## 6. Vessel layer

Recommended vessel IDs:

- `vessel-arm-vein`
- `vessel-right-heart`
- `vessel-pulmonary-artery`
- `vessel-pulmonary-vein`
- `vessel-left-heart`
- `vessel-arterial`
- `vessel-venous`
- `injection-site-left-arm`

Vessels should be centerline paths with no fill and a consistent stroke width.
React will apply concentration-dependent stroke colors and animated dash offset.

## 7. Delivery checks

Before replacing the placeholder map, verify:

```bash
xmllint --noout human-regional-map.svg
rg -o 'data-region="[^"]+"' human-regional-map.svg | sort
rg -o 'viewBox="[^"]+"' human-regional-map.svg
```

The second command should list all 24 region IDs exactly once at the interactive
group level. Also check that the SVG contains no absolute local file paths,
external raster URLs, scripts, or unsupported embedded HTML.

## 8. React integration

After delivery, the SVG paths replace the placeholder paths in
`app/organ-heatmap.tsx`. The ODE data, timeline, ranking, and color scale already
use the region IDs above, so no model recalculation is required merely to improve
the anatomy. The imported SVG must be inlined as React/SVG markup or transformed
into a React component; loading it only through an `<image>` tag would prevent
React from changing individual region fills.

The application currently uses Reactome's `Male body with organs`
(`R-ICO-013956`) twice: as a low-opacity reference layer and as same-origin SVG
fragment masks. Liver, lungs, kidneys, brain, and gastrointestinal anatomy are
therefore colored on their real source contours rather than on coarse blocks.
Custom aligned paths remain only where the source asset does not expose the 24
PBPK compartments separately. A future fully authored 24-region SVG can replace
the hybrid layer without changing ODE data, metric calculation, or selection
logic.
