class LabmethodsController < ApplicationController

  def index

  end

  def new

  end
  def create
      # find correct method to associate
      method_type= params["labmethod"]["methodtype"]
      selected_method = Methodtype.where(method_type: method_type)[0]

      # save the file
      uploaded_io = params["labmethod"]["pdf_file"]
      File.open(Rails.root.join('public','methods', method_type,
        uploaded_io.original_filename), 'wb') do |file|
          file.write(uploaded_io.read)
        end
      # grab params for labmethod
      title = params["labmethod"]["title"]
      pdf_file_name = uploaded_io.original_filename

      #create labmethod
      labmethod =Labmethod.new(title: title, pdf_url:pdf_file_name)
      labmethod.methodtype = selected_method

      labmethod.save
      redirect_to maintenance_index_path
  end

  def destroy

    labmethod = Labmethod.find(params[:id])
    labmethod.destroy
    redirect_to maintenance_index_path
  end
end
