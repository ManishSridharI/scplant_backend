import pandas as pd
import plotly.express as px
from django.http import JsonResponse
import numpy as np
from plotly.utils import PlotlyJSONEncoder
import scanpy as sc
import json

from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.JobAnnotateAndPlotModel import JobAnnotateAndPlotModel


@api_view(['POST'])
def DimPlot(request):
    try:
    # Retrieve method and annotation from the request with defaults
        method = request.data.get('method', 'UMAP').upper()
        annotation = request.data.get('annotation', 'Condition')
        job_annotate_and_plot_id = request.data.get('job_annotate_and_plot_id')
        
        job_annotate_and_plot_instance = JobAnnotateAndPlotModel.objects.get(id=job_annotate_and_plot_id)
        job_annotate_and_plot_file_output_instance = job_annotate_and_plot_instance.job_annotate_and_plot_file_output
        filepath = job_annotate_and_plot_file_output_instance.job_annotate_and_plot_output_with_celltype_file

        # Supported methods and corresponding columns in AnnData
        method_map = {
            'UMAP': 'X_umap',
            'T-SNE': 'X_tsne'
        }

        # Validate the method
        if method not in method_map:
            return JsonResponse({"error": f"Unsupported method: '{method}'. Supported methods are {list(method_map.keys())}"}, status=400)

        # Load AnnData object from file
        # filepath = "/home/sdkqc/Havel_portal/haval_backend/media/output_file/annotate_and_plot/output_with_celltype.h5ad"
        adata = sc.read_h5ad(filepath)

        # Extract coordinates dynamically based on the method
        coord_key = method_map[method]
        if coord_key not in adata.obsm:
            return JsonResponse({"error": f"Coordinates for '{method}' not found in AnnData."}, status=400)

        x_col, y_col = adata.obsm[coord_key][:, 0], adata.obsm[coord_key][:, 1]

        # Filter columns with <= 100 unique values for valid annotations
        valid_columns = [col for col in adata.obs.columns if adata.obs[col].nunique() <= 100]

        # Validate the requested annotation
        if annotation not in valid_columns:
            return JsonResponse({
                "error": f"Invalid annotation: '{annotation}'. Choose from {valid_columns}"
            }, status=400)

        # Prepare data for visualization
        data = pd.DataFrame({
            f"{method} 1": x_col,
            f"{method} 2": y_col,
            'Annotation': adata.obs[annotation].astype(str)
        })

        # Determine axis ranges
        x_min, x_max = data[f"{method} 1"].min(), data[f"{method} 1"].max()
        y_min, y_max = data[f"{method} 2"].min(), data[f"{method} 2"].max()

        # Create the scatter plot
        fig = px.scatter(
            data,
            x=f"{method} 1",
            y=f"{method} 2",
            color='Annotation',
            title=f"DimPlot {method} by {annotation}",
            labels={
                f"{method} 1": f"{method} 1",
                f"{method} 2": f"{method} 2",
                'Annotation': annotation
            },
            color_discrete_sequence=px.colors.qualitative.Light24
        )

        # Set all points visible by default with full opacity
        fig.update_traces(marker=dict(opacity=1))

        # Fix the axis ranges and tick marks
        fig.update_layout(
            margin=dict(r=210),
             modebar=dict(
        orientation='v',  # Horizontal toolbar on top
        # remove = ['Autoscale', 'zoomIn2d', 'zoomOut2d', 'lasso2d']
    ),
        xaxis=dict(
            range=[x_min, x_max],  # Lock the x-axis range
            showticklabels=True,  # Remove tick labels
        showgrid=False,        # Optionally remove the grid
        ticks="",              # Remove tick marks
    ),
         yaxis=dict(
        range=[y_min, y_max],  # Lock the y-axis range
        showticklabels=True,  # Remove tick labels
        showgrid=False,        # Optionally remove the grid
        ticks="",              # Remove tick marks
    ),
        legend=dict(
        x=1.1,  # Position to the right
        y=1,  # Align to the top
        traceorder="normal",
        orientation="v",  # Horizontal legend
        font=dict(size=12),  # Adjust legend font size
        bordercolor="Black",
        borderwidth=1,
        itemsizing="constant",
        itemclick="toggleothers",  # Clicking isolates the cluster
        itemdoubleclick="toggle"  # Double-click toggles visibility

    ),
        plot_bgcolor="white",  # Set the plot area background to white
        paper_bgcolor="white", 
  
        )

        data = [trace.to_plotly_json() for trace in fig.data]
        layout = fig.layout.to_plotly_json()
        # Manually handle any remaining ndarrays in traces
        for trace in data:
            if 'x' in trace and isinstance(trace['x'], (np.ndarray, pd.Series)):
                trace['x'] = trace['x'].tolist()
            if 'y' in trace and isinstance(trace['y'], (np.ndarray, pd.Series)):
                trace['y'] = trace['y'].tolist()



        # Return the serialized figure data and layout
        return JsonResponse({"data": data, "layout": layout, "columns": valid_columns})

    except Exception as e:
        return JsonResponse({"error happened": str(e)}, status=500)


@api_view(['POST'])
def DotPlot(request):
    try:
        # Retrieve annotation from the request with defaults
        annotation = request.data.get('annotation', 'scPlantAnnotate_celltype')
        marker  = request.data.get('marker', 'Top 5')
        job_annotate_and_plot_id = request.data.get('job_annotate_and_plot_id')
        
        job_annotate_and_plot_instance = JobAnnotateAndPlotModel.objects.get(id=job_annotate_and_plot_id)
        job_annotate_and_plot_file_output_instance = job_annotate_and_plot_instance.job_annotate_and_plot_file_output
        
        filepath = job_annotate_and_plot_file_output_instance.job_annotate_and_plot_output_with_celltype_file
        
        # Load AnnData object from file
        #filepath = "/home/sdkqc/Havel_portal/haval_backend/media/output_file/annotate_and_plot/output_with_celltype.h5ad"
        adata = sc.read_h5ad(filepath)


        # Validate the requested annotation
        if annotation not in adata.obs.columns:
            return JsonResponse({
                "error": f"Invalid annotation: '{annotation}'. Choose from {list(adata.obs.columns)}"
            }, status=400)
        

        valid_columns = [col for col in adata.obs.columns if adata.obs[col].nunique() <= 100]


        marker_filepath = {
                    "Top 5": job_annotate_and_plot_file_output_instance.job_annotate_and_plot_top5_markers_file,
                    "Top 10": job_annotate_and_plot_file_output_instance.job_annotate_and_plot_top10_markers_file,
                    "Top 25": job_annotate_and_plot_file_output_instance.job_annotate_and_plot_top25_markers_file
                }.get(marker, job_annotate_and_plot_file_output_instance.job_annotate_and_plot_top5_markers_file)  # Default to "Top 5" if marker is invalid

        # gettin the marker files
        # marker_filepath = {
        #     "Top 5" :"/home/sdkqc/Havel_portal/haval_backend/media/output_file/annotate_and_plot/top5_DEGs.xlsx",
        #     "Top 10" :"/home/sdkqc/Havel_portal/haval_backend/media/output_file/annotate_and_plot/top10_DEGs.xlsx",
        #     "Top 25" : "/home/sdkqc/Havel_portal/haval_backend/media/output_file/annotate_and_plot/top25_DEGs.xlsx"
        # }

        if marker not in marker_filepath:
            return JsonResponse({"error": f"Invalid marker: '{marker}'. Supported values are {list(marker_filepath.keys())}"}, status=400)

        marker_genes_filepath = marker_filepath[marker]
        marker_genes_df = pd.read_excel(marker_genes_filepath)
        marker_genes = marker_genes_df['gene'].tolist()

        # Extract data
        raw_adata = adata.raw.to_adata()  # Get raw data
        expression_data = raw_adata[:, marker_genes].X.toarray()
        expression_data = pd.DataFrame(expression_data, columns=marker_genes, index=adata.obs[annotation])

        # Compute average expression
        average_expression = expression_data.groupby(level=0, observed=False).mean()
        average_expression = average_expression.reset_index()
        average_expression = average_expression.melt(id_vars=annotation, var_name='Gene', value_name='Average Expression')

        # Compute fraction of cells expressing the gene
        expression_data_bool = expression_data.astype(bool)
        fraction_expression = expression_data_bool.groupby(level=0, observed=False).sum() / expression_data_bool.groupby(level=0, observed=False).count()
        fraction_expression = fraction_expression.reset_index()
        fraction_expression = fraction_expression.melt(id_vars=annotation, var_name='Gene', value_name='Fraction Expression')

        # Merge average expression and fraction data
        plot_data = pd.merge(average_expression, fraction_expression, on=[annotation, 'Gene'], how='inner')

        # Create the dot plot
        fig = px.scatter(
            plot_data,
            x='Gene',
            y=annotation,
            size='Fraction Expression',
            color='Average Expression',
            color_continuous_scale='YlOrBr',
            title='Dot Plot of Marker Genes',
        )

        # Update layout for better visualization
        fig.update_layout(
            xaxis_title='Genes',
            yaxis_title=annotation,
            coloraxis_colorbar=dict(title='Avg Expression'),
            plot_bgcolor='white',  # Remove background
            paper_bgcolor='white',
            margin=dict(b=60),
            xaxis=dict(tickangle=90),  # Rotate gene names for better readability
            yaxis=dict(autorange="reversed"),  # Reverse y-axis to match Scanpy's style
            legend=dict(
                x=1.1,  # Position legend to the right
                y=1,  # Align to the top
                traceorder="normal",
                font=dict(size=12),  # Adjust legend font size
                bordercolor="Black",
                borderwidth=1,
            ),
        )


        fig_json = json.loads(json.dumps(fig, cls=PlotlyJSONEncoder))

        data = fig_json["data"]
        layout = fig_json["layout"]

        # # Convert `fig.data` and `fig.layout` to JSON-serializable objects
        # data = [trace.to_plotly_json() for trace in fig.data]
        # layout = fig.layout.to_plotly_json()

        # Manually handle any remaining ndarrays in traces
        for trace in data:
            if 'x' in trace and isinstance(trace['x'], (np.ndarray, pd.Series)):
                trace['x'] = trace['x'].tolist()
            if 'y' in trace and isinstance(trace['y'], (np.ndarray, pd.Series)):
                trace['y'] = trace['y'].tolist()

        # Return the serialized figure data and layout
        return JsonResponse({"data": data, "layout": layout, 'columns': valid_columns})

    except Exception as e:
        return JsonResponse({"error happened": str(e)}, status=500)






